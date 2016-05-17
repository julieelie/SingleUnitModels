function pool = parpool(varargin)
%PARPOOL Create a parallel pool of workers on a cluster and return a pool object
%   pool = PARPOOL creates and returns a pool on the default cluster with its
%   NumWorkers in the range [1, preferredNumWorkers] for running parallel
%   language features (parfor, parfeval, parfevalOnAll, spmd, and distributed).
%   preferredNumWorkers is the value defined in your parallel preferences. Use
%   delete(pool) to shut down the parallel pool.
%
%   pool = PARPOOL(numWorkers) creates and returns a pool with the
%   specified number of workers. numWorkers can be a positive integer or a
%   range specified as a 2-element vector of integers. If numWorkers is a
%   range, the resulting pool has size as large as possible in the range
%   requested.
%
%   pool = PARPOOL(cluster) creates and returns a pool on the specified
%   cluster object.
%
%   pool = PARPOOL(profile) creates and returns a pool on the cluster
%   defined in the specified parallel profile.
%
%   pool = PARPOOL(profileOrCluster, numWorkers) creates and returns a pool
%   on the specified cluster or profile with the specified number of
%   workers.
%
%   pool = PARPOOL(..., P1, V1, ..., Pn, Vn) allows additional
%   parameter-value pairs that specify further properties of the created
%   pool. The accepted parameters are:
%
%   - 'AttachedFiles' - a string or a cell array of strings specifying
%   files to be attached to the pool at the time of its creation. Attach
%   files to the pool to allow them to be used by the workers.
%
%   - 'SpmdEnabled' - a logical value (true or false) that enables the use
%   of the SPMD language construct. By default 'SpmdEnabled' is set to
%   true. If SpmdEnabled is set to false then you will not be able to use
%   the SPMD language construct, however the parallel pool created should be
%   more robust against workers crashing.
%
%   - 'IdleTimeout' - an integer greater than zero that sets the time in
%   minutes after which the pool will automatically shut itself down if it
%   is idle. A pool is idle if it is not running code on the workers. By
%   default 'IdleTimeout' is the same as the value in your parallel preferences.
%
%   Examples:
%   % Create a pool which uses all the default settings.
%   myPool = parpool();
%
%   % Create a pool using the profile named 'MyProfile'.
%   myPool = parpool('MyProfile')
%
%   % Create a pool using the cluster myCluster.
%   % First create myCluster using 'parcluster'.
%   myCluster  = parcluster();
%   % Now create the pool.
%   myPool = parpool(myCluster);
%
%   % Create a pool with 2 workers.
%   myPool = parpool(2);
%
%   % Create a pool using default settings and attach the file 'myFunction.m'
%   % to it.
%   myPool = parpool('AttachedFiles', {'myFunction.m'});
%
%   % Create a pool using default settings and disable the use of SPMD.
%   myPool = parpool('SpmdEnabled', false);
%
%   % Create a pool using default settings and an IdleTimeout of 2 hours.
%   myPool = parpool('IdleTimeout', 120);
%
%   See also gcp, parfor, spmd, parfeval, parfevalOnAll,
%            parcluster, parallel.Cluster/parpool, batch,
%            parallel.Pool/addAttachedFiles,
%            parallel.Pool/updateAttachedFiles,
%            parallel.Pool/listAutoAttachedFiles,
%            parallel.Pool/delete.

%   Copyright 2013-2014 The MathWorks, Inc.

% Since we do not yet have a pool and therefore no Session or
% SessionInfo we need to create a NULL_SESSION_INFO in order to notify
% Java listeners in the event of any errors.
sessionInfo = com.mathworks.toolbox.distcomp.pmode.SessionInfo.NULL_SESSION_INFO;
% This cleanupObj will notify listeners of the termination of this function
% by the user.
cleanupObj = parallel.internal.general.DisarmableOncleanup(...
    @() parallel.internal.pool.notifyAbortByUser(sessionInfo));

try
    pool = parallel.internal.pool.doParpool(varargin{:});
catch err
    % The following nested try-catch block is necessary in order that the
    % error passed to the SessionFactory notifier includes parpool in the stack trace
    try
        throw(err)
    catch parpoolErr
        if ~strcmp(parpoolErr.identifier, 'parallel:convenience:ConnectionOpen')
            parallel.internal.pool.notifySessionFailedToStart(sessionInfo, parpoolErr);
        end
        % Since a more specific notification can be made, we can disarm the
        % original cleanupObj
        cleanupObj.disarm();        
    end    
    throw(err);
end
% Now completed pool creation, the cleanupObj can be disarmed
cleanupObj.disarm();
end
