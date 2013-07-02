% SQL_fdep
% Finds dependencies between pieces of code in OperationCode using fdep
% Also finds Matlab toolboxes and periphery toolboxes and labels all time series accordingly
% Writes the results to mySQL database (CodeSource, CodeSourceLink)
% Ben Fulcher July 2013

getfromwhere = 'CodeTable'; % 'opsdir'

% Open database connection:
[dbc, dbname] = SQL_opendatabase;

%% RETRIEVE Operations to analyze
switch getfromwhere
case 'opsdir' % Get operation list from Operations directory
    % Directory containing all the operations
    theopsat = '/Users/ben/Dropbox/Work/Time Series Documentation/HCTSpackage/Operations';

    %% Retrieve all .m files from the directory
    a = dir(theopsat); % read the directory
    filenames = {a.name}; % cell of filenames of operations
    nfiles = length(filenames);
    exts = cell(nfiles,1);
    for i = 1:nfiles
        [~, exts{i}] = strtok(filenames{i},'.');
    end
    wewant = strcmp('.m',exts); % only retrieve .m files
    filenames = filenames(wewant);
    
case 'CodeTable' % Get operation list from OperationCode table in database
    SelectString = 'SELECT c_id, CodeName FROM OperationCode WHERE fdepDone = 0'; % WHERE//later you can only calculate new ones
    [qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);
    if ~isempty(emsg), beep; keyboard; end
    cids_anal = vertcat(qrc{:,1}); % the c_ids analyzed in this run of the script
    filenames = qrc(:,2);
    
end

nfiles = length(filenames); % number of files found

%% 1. Save screen output to text file
% fn = 'Operation_fdep.txt'; % filename
% fid = fopen(fn,'w');
% for i = 1:nfiles
%     fprintf(fid,'----------_---------%s----------_---------\n',filenames{i})
%     T = evalc('fdep(filenames{i});'); % get text output
%     fprintf(fid,'%s\n\n',T);
% end

%% Perform more subtle analysis to store back to database
fid = 1; % write results to screen

toolboxes = cell(nfiles,1); % save toolboxes used by each method
modules = cell(nfiles,1); % save modules used 1-level down, by each method
downstream = cell(nfiles,1); % save all downstream modules used by each method
dashes = '----------------------';
fprintf(1,'Running fdep to find functional dependencies between %u pieces of Operation code\n',nfiles)
fdepworked = logical(ones(nfiles,1));
for i = 1:nfiles
    fprintf(fid,'%s%s [%u/%u]%s\n',dashes(1:end-ceil(length(filenames{i})/2)),filenames{i},i,nfiles,dashes(1:end-floor(length(filenames{i})/2)));
    try
        [T, p] = evalc('fdep(filenames{i})'); % get output
        if isempty(p)
            fprintf(1,'fdep provided no output for %s\n',filenames{i})
            fdepworked(i) = 0;
            continue
        end
    catch
        fdepworked(i) = 0;
        fprintf(1,'Error running fdep on %s\n',filenames{i});
    end
        
    toolboxnow = p.toolbox(cellfun(@(x)any(x==1),p.modbix));
                    % toolboxes used DIRECTLY by the input function
                    % NB: periphery files that use the toolbox are not included
    if ~isempty(toolboxnow)
        components = cellfun(@(x)regexp(x,'\[(\w*)\]','match'),toolboxnow,'UniformOutput',0);
        components = cellfun(@(x)x{1},components,'UniformOutput',0); % change cell of cells to cell of strings
        components = cellfun(@(x)[x(2:end-1),' toolbox'],components,'UniformOutput',0); % remove the surrounding square brackets '[ ]'
        toolboxes{i} = components; % cell of strings containing matlab toolboxes used DIRECTLY in this code 
                                                % (i.e., doesn't include subsidiary codes)
    else
        toolboxes{i} = '';
    end
    modules{i} = p.module(p.depth(:,2) == 2); % only store one-level down dependencies,
                                                % subsequent levels stored later
    downstream{i} = p.module(p.depth(:,2) > 1); % store all downstream user-defined functions, to add to Code table
end

% Find auxiliary files and add them to the OperationCode Table in the database
SelectString = 'SELECT CodeName FROM OperationCode';
codenames_db = mysql_dbquery(dbc,SelectString); % all codenames in the database

allmodules = sort(unique(vertcat(downstream{:}))); % all user functions called by code
if ~isempty(allmodules)
    isnew = ~ismember(allmodules,codenames_db);
    nnew = sum(isnew);
else
    nnew = 0;
end
if nnew > 0 % some new dependent subsidiary operation files need to be added to the database
    % display the new ones
    fprintf(1,'The following %u subsidiary files are called by existing files in the database. I need to add them to %s\n',nnew,dbname);
    fisnew = find(isnew); % indicies of new OperationCodes
    for i = 1:nnew
        fprintf(1,'%s\n',allmodules{fisnew(i)})
    end
    input('ok??')
    esc = @sqlescapestring; % inline function to add escape strings to format mySQL queries
    for i = 1:nnew
        InsertString = sprintf('INSERT INTO OperationCode (CodeName, IsSubsidiary, fdepDone) VALUES (''%s'',1,0)',esc(allmodules{fisnew(i)}));
        [~,emsg] = mysql_dbexecute(dbc,InsertString);
        if ~isempty(emsg), fprintf(1,'Error with %s\n',allmodules{fisnew(i)}); keyboard; end
    end
    % Ok, so now we have all of the operations used by all of the other operations(!)
    % (may need to rerun to get all the way down the hierarchy??)
    %% REALLY NEED THIS TO BE A WHILE LOOP FOR WHILE nadded > 0
else
    fprintf(1,'No new files to add :))\n')
end


% Fill the CodeLinks table -- should link every piece of code with immediate links referenced in that code
SelectString = 'SELECT c_id, CodeName FROM OperationCode';
[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);
cids = vertcat(qrc{:,1});
codenames_db = qrc(:,2);

addcell = {};
for i = 1:nfiles
    % look at the modules and add elements
    if ~isempty(modules{i}) % relies on external code
        cid_i = cids(strcmp(codenames_db,filenames{i})); % the c_id of file i
        if isempty(cid_i); fprintf(1,'empty at %s\n',filenames{i});
            continue; % (this operation is in the directory but never referenced)
        elseif length(cid_i) > 1
            fprintf(1,'%s referenced %u times??!\n',filenames{i},length(cid_i));
            keyboard
        end
        for j = 1:length(modules{i})
            cid_j = cids(strcmp(codenames_db,modules{i}{j})); % the c_id of file i
            addcell{end+1} = sprintf('(%u,%u)',cid_i,cid_j); % yes it grows within a loop so sue me
        end
    end
end
SQL_add_chunked(dbc,'INSERT INTO OperationCodeLinks (c_id_up, c_id_down) VALUES',addcell); % add them all in chunks
fprintf(1,'Added %u new links between operation code files to %s\n',length(addcell),dbname)


%% ASSIGN TOOLBOXES/CodeSources
% CodeSource table
% CodeSourceLink table
% First we can add the Matlab toolboxes found by fdep
distinct_tb = unique(vertcat(toolboxes{:})); % all toolboxes used by the current set of operations
if ~isempty(distinct_tb) % some toolboxes registered -- check if they're new
    % find new ones:
    codesourcenames = mysql_dbquery(dbc,'SELECT Name FROM CodeSource');
    isduplicate = ismember(distinct_tb,codesourcenames);
    if any(~isduplicate) % some new toolbox sources to add to the database
        fprintf(1,'%u new toolbox sources to add to CodeSource...',sum(~isduplicate))
        addcell = {};
        for i = 1:length(distinct_tb)
            if ~isduplicate(i) % add this to the CodeSource table
                addcell{end+1} = sprintf('(''%s'')',distinct_tb{i});
            end
        end
        SQL_add_chunked(dbc,'INSERT INTO CodeSource (Name) VALUES',addcell);
        fprintf(1,' added.\n')
    end
end

% Add CodeSource links
% for each of the nfiles, we need to add links to toolboxes in CodeSourceLinks
qrc = mysql_dbquery(dbc,'SELECT tb_id, Name FROM CodeSource');
tbids = vertcat(qrc{:,1}); % all toolbox ids in CodeSource table
tbnames_db = qrc(:,2);
addcell = {};
for i = 1:nfiles
    if ~isempty(toolboxes{i}) % this file used toolboxes
        cid_i = cids(strcmp(codenames_db,filenames{i})); % the c_id of file i
        if isempty(cid_i); fprintf(1,'~~~empty at %s\n',filenames{i});
            continue; % (this operation is in the directory but never referenced)
        elseif length(cid_i) > 1
            fprintf(1,'%s referenced %u times??!\n',filenames{i},length(cid_i));
            keyboard
        end
        for j = 1:length(toolboxes{i})
            tbid_j = tbids(strcmp(tbnames_db,toolboxes{i}{j})); % the tb_id of this toolbox name
            addcell{end+1} = sprintf('(%u,%u)',cid_i,tbid_j); % yes it grows within a loop so sue me
        end
    end
end
SQL_add_chunked(dbc,'INSERT INTO CodeSourceLink (c_id, tb_id) VALUES',addcell); % add them all in chunks
fprintf(1,'%u new links added between operation code files and Matlab toolboxes to CodeSourceLink table of %s\n',length(addcell),dbname)

%% Now we need to assign rules for other toolboxes
fprintf(1,'Now assigning links to different toolbox sources using directory information\n')
% NB: Only does this for the initial set of operations, not any subsidiaries added since
addsource = cell(nfiles,1);
for i = 1:nfiles
    whereami = which(filenames{i});
    hier_dir = regexp(whereami,'/','split'); % set of directories leading to HCTS package
    here = find(strcmp(hier_dir,'HCTSpackage')); % index of HCTSpackage directory
    if sum(here) == 0 % HCTSpackage doesn't appear in this path??
        fprintf(1,'%s isn''t in HCTSpackage directory??\n');
        addsource{i} = 'external';
    else
        nextloc = hier_dir{here+1};
        switch nextloc
        case 'Operations'
            % Code files in Operations directory are assigned to a source 'hctsa'
            addsource{i} = 'hctsa';
            fprintf(1,'%s is in the Operations directory of HCTSA\n',filenames{i})
        case 'PeripheryFunction'
            % Code files in PeripheryFunctions are assigned to a source 'hctsa-periphery'
            addsource{i} = 'hctsa-periphery';
            fprintf(1,'%s is in the PeripheryFunctions directory of HCTSA\n',filenames{i})
        case 'Toolboxes'
            % Code files in Toolboxes/xxx/ are assigned to a source 'xxx'
            thetoolbox = hier_dir{here+2};
            addsource{i} = thetoolbox; % the name of the directory
            fprintf(1,'%s is in the Toolbox directory %s of HCTSA\n',filenames{i},thetoolbox)
        otherwise % unknown place
            addsource{i} = ''; % we don't know where it is
        end
    end
end

% So now we add unique addsources to the CodeSource table that aren't duplicates
alldirsources = sort(unique(addsource)); % tbnames_db
isduplicate = ismember(alldirsources,tbnames_db);
if any(~isduplicate) % some new ones, add them
    for i = 1:length(alldirsources)
        if ~isduplicate(i) % a new one, add it
            [~,emsg] = mysql_dbexecute(dbc,sprintf('INSERT INTO CodeSource (Name) VALUES (''%s'')',alldirsources{i}));
            fprintf(1,'Added %s to the CodeSource table in %s\n',alldirsources{i},dbname)
        end
    end
    % has changed since the new insert -- will be needed for the linking performed below
    qrc = mysql_dbquery(dbc,'SELECT tb_id, Name FROM CodeSource');
    tbids = vertcat(qrc{:,1}); % all toolbox ids in CodeSource table
    tbnames_db = qrc(:,2);
end
% Now add all the links, and we'll (finally) be done!!

addcell = {};
for i = 1:nfiles
    if ~isempty(addsource{i}) % this file needs a directory-based codesource added
        cid_i = cids(strcmp(codenames_db,filenames{i})); % the c_id of file i
        if isempty(cid_i); fprintf(1,'empty at %s\n',filenames{i});
            continue; % (this operation is in the directory but never referenced)
        elseif length(cid_i) > 1
            fprintf(1,'%s referenced %u times??!\n',filenames{i},length(cid_i));
            keyboard
        end
        tbid_i = tbids(strcmp(tbnames_db,addsource{i})); % the tb_id of this source name
        addcell{end+1} = sprintf('(%u,%u)',cid_i,tbid_i); % yes it grows within a loop so sue me
    end
end
SQL_add_chunked(dbc,'INSERT INTO CodeSourceLink (c_id, tb_id) VALUES',addcell); % add them all in chunks
fprintf(1,'Added %u new links between operation code files and directory-derived sources to CodeSourceLink table of %s\n',length(addcell),dbname)


% Note the files analyzed using the fdepDone tag in OperationCode table
UpdateString = sprintf('UPDATE OperationCode SET fdepDone = 1 WHERE c_id IN (%s)',bencat(cids_anal(fdepworked))); % the c_ids for which fdep worked
[~,emsg] = mysql_dbexecute(dbc,UpdateString);

% Update the Ncodes count in CodeSource by counting references in the linking table
for i = 1:length(tbids)
    UpdateString = sprintf('UPDATE CodeSource SET Ncodes = (SELECT COUNT(*) FROM CodeSourceLink WHERE tb_id = %u) WHERE tb_id = %u',tbids(i),tbids(i));
    [~,emsg] = mysql_dbexecute(dbc,UpdateString);
end

fprintf(1,'%u pieces of operation code were analyzed and updated; this is now reflected in the database\n',length(cids_anal))

% Check whether there are pieces of code with fdep = 0
qrc = mysql_dbquery(dbc,'SELECT COUNT(*) FROM OperationCode WHERE fdepDone = 0');
fprintf(1,'Ok done. There are %u pieces of code still in OperationCode with unknown file dependencies. You could re-run this...\n',qrc{1})

SQL_closedatabase(dbc)