# create table of stage names 

drop table if exists _stats_stagenames;
create temporary table _stats_stagenames as
    select stageid, stagename  
    from Logs
    group by stagename
    order by stageid
    ;    
    
# stats per visit

drop table if exists _stats_visit;
create temporary table _stats_visit as 
    select loopnum, workerid, sum(duration) as dur
    from Durations where name like "harness.pipeline.visit.stage" 
    group by loopnum, workerid
    ;
        
drop table if exists _stats_visit_app;
create temporary table _stats_visit_app as 
    select loopnum, workerid, sum(duration) as dur
    from Durations where name like "harness.slice.visit.stage.process" 
    group by loopnum, workerid
    ;
    
select count( * ) as N,
        avg( p.dur )/1000000000 as average_app,
        std( p.dur )/1000000000 as stddev_app,
        avg( v.dur - p.dur )/1000000000 as average_pipe,
        std( v.dur - p.dur )/1000000000 as stddev_pipe,
        avg( v.dur )/1000000000 as average_visit,
        std( v.dur )/1000000000 as stddev_visit
    from _stats_visit_app as p, _stats_visit as v
    where p.loopnum = v.loopnum
        and p.workerid = v.workerid
    ;

# create the temporary table with stage stats;
# each row in this table represents a single stage
# in a single visit.

# Subtlety 1: there will be one extra "preprocess" per
# workerid for setting up a worker when the data has
# been exhausted. There is no "postprocess" for those
# cases. That is why I'm starting this table based on
# available "postprocess" data, which will result in the
# superfluous "preprocess" data being ignored.

# Subtlety 2: couldn't find a way to tell the select
# statement to create the placeholder columns pre_dur
# and app_dur as bigint(20), so these columns are 
# explicitly altered after the table is created.

# Subtlety 3: rather than doing big multidimensional
# joins on loopnum, workerid, and stageid simultaneously,
# it makes more sense to create a column I call "hashkey"
# that contains all that information and then gets indexed
# for speed; joins are then done on that indexed field.

# pre_dur = "harness.pipeline.visit.stage.preprocess"
# app_dur = "harness.slice.visit.stage.process"
# post_dur = "harness.pipeline.visit.stage.postprocess"
# stage_dur = "harness.pipeline.visit.stage"

drop table if exists _stats_stage;
create temporary table _stats_stage as
    select loopnum, workerid, stageid,
        duration as post_dur, -1 as pre_dur, -1 as app_dur,
        -1 as stage_dur,
        concat( workerid, "/", loopnum, "/", stageid ) as hashkey
    from Durations
    where name like "harness.pipeline.visit.stage.postprocess" 
    group by loopnum, workerid, stageid
    ;
alter table _stats_stage modify pre_dur bigint(20);
alter table _stats_stage modify app_dur bigint(20);
alter table _stats_stage modify stage_dur bigint(20);
alter table _stats_stage add index ( hashkey );

# get the information on app time so it can be merged
# into _stats_stage

drop table if exists _stats_app_by_stage;
create temporary table _stats_app_by_stage as
    select loopnum, workerid, stageid, duration as dur,
        concat( workerid, "/", loopnum, "/", stageid ) as hashkey
    from Durations
    where name like "harness.slice.visit.stage.process"
    group by loopnum, workerid, stageid
    ;
alter table _stats_app_by_stage add index ( hashkey );

# merge app timing information into _stats_stage    
    
update _stats_stage t1
    join _stats_app_by_stage app 
        on app.hashkey = t1.hashkey
    set t1.app_dur = app.dur
    ;
drop table _stats_app_by_stage;

# get the information on slice time so it can be merged
# into _stats_stage

drop table if exists _stats_slice_by_stage;
create temporary table _stats_slice_by_stage as
    select loopnum, workerid, stageid, duration as dur,
        concat( workerid, "/", loopnum, "/", stageid ) as hashkey
    from Durations
    where name like "harness.pipeline.visit.stage"
    group by loopnum, workerid, stageid
    ;
alter table _stats_slice_by_stage add index ( hashkey );

# merge app timing information into _stats_stage    
    
update _stats_stage t1
    join _stats_slice_by_stage slice 
        on slice.hashkey = t1.hashkey
    set t1.stage_dur = slice.dur
    ;
drop table _stats_slice_by_stage;

# get the information on pre-process time
# so it can be merged into _stats_stage

drop table if exists _stats_preproc_by_stage;
create temporary table _stats_preproc_by_stage as
    select loopnum, workerid, sum(duration) as dur, stageid,
        concat( workerid, "/", loopnum, "/", stageid ) as hashkey
    from Durations where name like "harness.pipeline.visit.stage.preprocess" 
        and stageid > 0
    group by loopnum, workerid, stageid
    ;    
alter table _stats_preproc_by_stage add index ( hashkey );

# merge preprocess timing information into _stats_stage

update _stats_stage t1
    join _stats_preproc_by_stage pre 
        on pre.hashkey = t1.hashkey
    set t1.pre_dur = pre.dur
    ;
drop table _stats_preproc_by_stage;

# _stats_stage now fully loaded

select count( * ) as N,
    avg( stage_dur - app_dur ) / 1000000000 as avg_pipe,
    std( stage_dur - app_dur ) / 1000000000 as std_pipe
    from _stats_stage
    ;

select count( * ) as N,
    avg( stage_dur - app_dur ) / 1000000000 as avg_pipe,
    std( stage_dur - app_dur ) / 1000000000 as std_pipe
    from _stats_stage
    where stageid > 1 and stageid < 290
    ;    

# get all stage stats

select count( * ) as N,
    n.stageid as stageid,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app,
    avg( s.stage_dur - s.app_dur ) / 1000000000 as avg_pipe,
    std( s.stage_dur - s.app_dur ) / 1000000000 as std_pipe
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid
    group by stage
    order by stageid
    ;
   
# grouping some stages together

select count( * ) as N,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid and
        n.stagename like "isrInput%" and
        n.stagename not like "isrInputRaw"
    order by n.stageid
    ;   
   
select count( * ) as N,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid and
        n.stagename like "isrOutput%"
    order by n.stageid
    ;   
   

select count( * ) as N,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid and
        n.stagename like "isrSaturation%"
    order by n.stageid
    ;

select count( * ) as N,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid and
        n.stagename like "isrOverscan%"
    order by n.stageid
    ;
    
select count( * ) as N,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid and
        n.stagename like "isrBias%"
    order by n.stageid
    ;
    
select count( * ) as N,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid and
        n.stagename like "isrVariance%"
    order by n.stageid
    ;
     
select count( * ) as N,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid and
        n.stagename like "isrDark%"
    order by n.stageid
    ;
    
select count( * ) as N,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid and
        n.stagename like "isrFlat%"
    order by n.stageid
    ;
   
select count( * ) as N,
    n.stagename as stage,
    avg( s.app_dur ) / 1000000000 as avg_app,
    std( s.app_dur ) / 1000000000 as std_app
    from _stats_stage as s,
        _stats_stagenames as n
    where n.stageid = s.stageid and
        n.stagename like "isrSdqaAmp%"
    order by n.stageid
    ;
   
    