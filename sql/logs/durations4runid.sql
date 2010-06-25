-- clean up from previous run
drop view startloop;
drop view runtbl;
set @rid = 'DG0149';


-- create the table
select min(micros) from logger where runid=@rid;
create table DG0149 as select *, (micros-1203371662663923)/1000000.0 
       as secs from logger where runid=@rid;
-- 187463 rows affected, 65535 warnings (6.08 sec)
alter table DG0149 add column loopnum int(5), add column stage int(3);
create view runtbl as select * from DG0149;

-- some indexes will help
create index nodeidx on DG0149 (node);
create index sliceidx on DG0149 (sliceid);

-- fill in the node column
update runtbl set 
       node=convert(substr(hostId from 5 for locate('.',hostId)-5), signed);
-- 187463 rows affected (14.32 sec)

create view startloop as select * from runtbl 
                         where comment like '%ing stage loop number%';

call fillallnodes();  -- (52 min 52.77 sec)
call fixendloop();    -- (9 min 16.46 sec)

-- calculating durations
-- loop:
insert into durations 
       select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
              t1.loopnum,t1.stage,'master loop' as name,t1.id,
              t2.secs-t1.secs as duration 
       from runtbl t1, runtbl t2 
       where t1.comment like 'starting stage loop%' and 
             t2.comment like 'starting stage loop%' and 
             t1.sliceid=-1 and t2.sliceid=-1 and t1.node=t2.node and 
             t1.loopnum+1=t2.loopnum 
       order by t1.node,t1.loopnum,t1.sliceid;
-- 186 rows affected (28.40 sec)

-- stage:
insert into durations 
       select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
              t1.loopnum,t1.stage,'stage' as name,t1.id,
              t2.secs-t1.secs as duration
       from runtbl t1, runtbl t2
       where t1.node=t2.node and t1.sliceid=-1 and t2.sliceid=-1 and 
             t1.comment like 'iStage%' and t2.comment like 'iStage%' and 
             t1.stage+1=t2.stage and t1.loopnum=t2.loopnum 
       order by t1.node,t1.loopnum,t1.stage;
-- 1178 rows affected (3 min 10.16 sec)
insert into durations 
       select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
              t1.loopnum,t1.stage,'stage' as name,t1.id,
              t2.secs-t1.secs as duration
       from runtbl t1, runtbl t2
       where t1.node=t2.node and t1.sliceid=-1 and t2.sliceid=-1 and 
             t1.comment like 'iStage%' and 
             t2.comment like 'starting stage loop%' and 
             t1.stage=10 and t1.node=6 and t2.node=6 and t1.loopnum+1=t2.loopnum 
       order by t1.node,t1.loopnum,t1.stage;
-- 62 rows affected (0.63 sec)
insert into durations 
       select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
              t1.loopnum,t1.stage,'stage' as name,t1.id,
              t2.secs-t1.secs as duration
       from runtbl t1, runtbl t2
       where t1.node=t2.node and t1.sliceid=-1 and t2.sliceid=-1 and 
             t1.comment like 'iStage%' and 
             t2.comment like 'starting stage loop%' and 
             t1.stage=4 and t1.node=9 and t2.node=9 and t1.loopnum+1=t2.loopnum 
       order by t1.node,t1.loopnum,t1.stage;
-- 62 rows affected (0.20 sec)
insert into durations 
       select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
              t1.loopnum,t1.stage,'stage' as name,t1.id,
              t2.secs-t1.secs as duration
       from runtbl t1, runtbl t2
       where t1.node=t2.node and t1.sliceid=-1 and t2.sliceid=-1 and 
             t1.comment like 'iStage%' and 
             t2.comment like 'starting stage loop%' and 
             t1.stage=8 and t1.node=10 and t2.node=10 and 
             t1.loopnum+1=t2.loopnum 
       order by t1.node,t1.loopnum,t1.stage;
-- 62 rows affected (0.26 sec)

-- process():
insert into durations 
        select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
               t1.loopnum,t1.stage,'stage process' as name,t1.id,
               t2.secs-t1.secs as duration 
        from runtbl t1, runtbl t2 
        where t1.comment like '%Starting process%' and 
              t2.comment like '%Ending process%' and 
              t1.sliceid=t2.sliceid and t1.node=t2.node and 
              t1.loopnum=t2.loopnum and t1.stage=t2.stage 
        order by t1.node,t1.loopnum,t1.sliceid;
-- 24056 rows affected (34 min 5.58 sec)

-- preprocess():
insert into durations 
        select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
               t1.loopnum,t1.stage,'stage preprocess' as name,t1.id,
               t2.secs-t1.secs as duration 
        from runtbl t1, runtbl t2 
        where t1.comment like '%Starting preprocess%' and 
              t2.comment like '%Ending preprocess%' and 
              t1.sliceid=-1 and t2.sliceid=-1 and t1.node=t2.node and 
              t1.loopnum=t2.loopnum and t1.stage=t2.stage 
        order by t1.node,t1.loopnum,t1.sliceid;
-- 1364 rows affected (3 min 38.67 sec)

-- postprocess():
insert into durations 
        select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
               t1.loopnum,t1.stage,'stage postprocess' as name,t1.id,
               t2.secs-t1.secs as duration 
        from runtbl t1, runtbl t2 
        where t1.comment like '%Starting postprocess%' and 
              t2.comment like '%Ending postprocess%' and 
              t1.sliceid=-1 and t2.sliceid=-1 and t1.node=t2.node and 
              t1.loopnum=t2.loopnum and t1.stage=t2.stage 
        order by t1.node,t1.loopnum,t1.sliceid;
-- 1364 rows affected (3 min 37.57 sec)

-- process wait:
insert into durations 
        select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
               t1.loopnum,t1.stage,'process sync wait' as name,t1.id,
               t2.secs-t1.secs as duration 
        from runtbl t1, runtbl t2 
        where t1.comment like '%Getting post-process signal%' and 
              t2.comment like '%Starting postprocess%' and 
              t1.sliceid!=-1 and t2.sliceid=-1 and t2.node<9 and t1.node<9 and
              t1.loopnum=t2.loopnum and t1.stage=t2.stage 
        order by t1.node,t1.loopnum,t1.sliceid;
-- 22320 rows affected (1 min 23.45 sec)
insert into durations 
        select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
               t1.loopnum,t1.stage,'process sync wait' as name,t1.id,
               t2.secs-t1.secs as duration 
        from runtbl t1, runtbl t2 
        where t1.comment like '%Getting post-process signal%' and 
              t2.comment like '%Starting postprocess%' and 
              t1.sliceid!=-1 and t2.sliceid=-1 and t2.node=9 and t1.node=9 and 
              t1.loopnum=t2.loopnum and t1.stage=t2.stage 
        order by t1.node,t1.loopnum,t1.sliceid;
-- 1240 rows affected (49.66 sec)
insert into durations 
        select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
               t1.loopnum,t1.stage,'process sync wait' as name,t1.id,
               t2.secs-t1.secs as duration 
        from runtbl t1, runtbl t2 
        where t1.comment like '%Getting post-process signal%' and 
              t2.comment like '%Starting postprocess%' and 
              t1.sliceid!=-1 and t2.sliceid=-1 and t2.node=10 and t1.node=10 and
              t1.loopnum=t2.loopnum and t1.stage=t2.stage 
        order by t1.node,t1.loopnum,t1.sliceid;
-- 496 rows affected (1 min 2.51 sec)

-- event wait:
insert into durations
       select t1.runid,t1.sliceid,t1.date,t1.micros,t1.node,t1.visitId,t1.secs,
               t1.loopnum,t1.stage,'event wait' as name,t1.id,
               t2.secs-t1.secs as duration 
       from runtbl t1, runtbl t2 
       where t1.comment like 'waiting on receive%' and 
             t2.comment like 'received event; sending%' and 
             t1.sliceid=-1 and t2.sliceid=-1 and t2.node=t1.node and
             t1.loopnum=t2.loopnum and t1.stage=t2.stage 
       order by t1.node,t1.loopnum,t1.sliceid,t1.stage;
-- 248 rows affected (39.39 sec)

-- loop event wait:
-- create view loopeventwait as
--        select runid,node,loopnum,sum(duration) as sum from durations 
--               where name='event wait' 
--               group by runid,node,loopnum;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'loop event wait' as name,d.id,s.sum as duration
       from durations d, loopeventwait s
       where d.stage=1 and d.loopnum=s.loopnum and d.node=s.node and 
             d.name='event wait' and d.runid=@rid and s.runid=@rid
       order by d.node,d.loopnum;
-- 186 rows affected (0.20 sec)

-- max process():
-- create view maxproctime as 
--        select runid,loopnum,stage,max(duration) as maxdur
--        from durations 
--        where name='stage process' and node!=9 and node!=10 
--        group by runid,loopnum,stage;
-- create view maxproc as 
--        select m.runid,m.loopnum,m.stage,m.maxdur,max(d.id) as maxid 
--        from maxproctime m, durations d 
--        where d.name='stage process' and m.runid=d.runid and 
--              m.loopnum=d.loopnum and m.stage=d.stage and 
--              m.maxdur=d.duration and d.node!=9 and d.node!=10
--        group by m.runid,m.loopnum,m.stage,m.maxdur;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'max stage process' as name,d.id,d.duration 
       from durations d, maxproc m 
       where d.loopnum=m.loopnum and d.stage=m.stage and 
             d.id=m.maxid and d.name='stage process' and 
             m.runid=@rid and d.runid=@rid
       order by loopnum,stage;
-- 620 rows affected (9.79 sec)
-- create view maxproc9time as 
--        select runid,loopnum,stage,max(duration) as maxdur
--        from durations 
--        where name='stage process' and node=9
--        group by runid,loopnum,stage;
-- create view maxproc9 as 
--        select m.runid,m.loopnum,m.stage,m.maxdur,max(d.id) as maxid 
--        from maxproc9time m, durations d 
--        where d.name='stage process' and m.runid=d.runid and d.node=9 and
--              m.loopnum=d.loopnum and m.stage=d.stage and m.maxdur=d.duration 
--        group by m.runid,m.loopnum,m.stage,m.maxdur;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'max stage process' as name,d.id,d.duration 
       from durations d, maxproc9 m 
       where d.loopnum=m.loopnum and d.stage=m.stage and d.node=9 and 
             d.id=m.maxid and d.name='stage process' and 
             m.runid=@rid and d.runid=@rid
       order by loopnum,stage;
-- 248 rows affected (0.43 sec)
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'max stage process' as name,d.id,d.duration 
       from durations d 
       where d.name='stage process' and d.node=10 and d.runid=@rid
       order by loopnum,stage;
-- 496 rows affected (0.09 sec)

--- app code:
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'app' as name,d.id,
              d.duration+s.duration+p.duration as duration
       from durations d, durations s, durations p 
       where d.name='stage preprocess' and  s.name='stage postprocess' and 
             p.name='max stage process' and 
             d.loopnum=s.loopnum and s.loopnum=p.loopnum and 
             d.stage=s.stage and s.stage=p.stage and 
             d.node=s.node and s.node<9 and p.node<9 and 
             d.runid=@rid and d.runid=s.runid and d.runid=p.runid
       order by d.loopnum,d.stage,d.node;
-- 620 rows affected (2.45 sec)
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'app' as name,d.id,
              d.duration+s.duration+p.duration as duration
       from durations d, durations s, durations p 
       where d.name='stage preprocess' and  s.name='stage postprocess' and 
             p.name='max stage process' and 
             d.loopnum=s.loopnum and s.loopnum=p.loopnum and 
             d.stage=s.stage and s.stage=p.stage and 
             d.node=s.node and s.node=p.node and s.node=9 and 
             d.runid=@rid and d.runid=s.runid and d.runid=p.runid
       order by d.loopnum,d.stage,d.node;
-- 248 rows affected (0.64 sec)
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'app' as name,d.id,
              d.duration+s.duration+p.duration as duration
       from durations d, durations s, durations p 
       where d.name='stage preprocess' and  s.name='stage postprocess' and 
             p.name='max stage process' and 
             d.loopnum=s.loopnum and s.loopnum=p.loopnum and 
             d.stage=s.stage and s.stage=p.stage and 
             d.node=s.node and s.node=p.node and s.node=10 and 
             d.runid=@rid and d.runid=s.runid and d.runid=p.runid
       order by d.loopnum,d.stage,d.node;
-- 496 rows affected (1.00 sec)

--- app code over entire loop:
-- create view apploop as 
--        select runid,node,loopnum,sum(duration) as sum from durations 
--        where name='app' group by runid,node,loopnum;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'app loop' as name,d.id,s.sum as duration
       from durations d, apploop s
       where d.stage=1 and d.loopnum=s.loopnum and d.node=s.node and 
             d.name='app' and d.runid=@rid and d.runid=s.runid
       order by d.node,d.loopnum;
-- 186 rows affected (0.20 sec)

-- mw overhead per stage:
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'stage overhead' as name,d.id,
              d.duration-s.duration-w.duration as duration
       from durations d, durations s, durations w
       where d.name='stage' and  s.name='app' and w.name='event wait' and 
             d.loopnum=s.loopnum and s.loopnum=w.loopnum and 
             d.stage=1 and d.stage=s.stage and s.stage=w.stage and 
             d.node=s.node and s.node=w.node and 
             d.runid=@rid and d.runid=s.runid and d.runid=w.runid
       order by d.node,d.loopnum,d.stage;
-- 186 rows affected (0.32 sec)
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'stage overhead' as name,d.id,
              d.duration-s.duration-w.duration as duration
       from durations d, durations s, durations w
       where d.name='stage' and  s.name='app' and w.name='event wait' and 
             d.loopnum=s.loopnum and s.loopnum=w.loopnum and 
             d.stage=5 and d.stage=s.stage and s.stage=w.stage and 
             d.node=10 and d.node=s.node and s.node=w.node and 
             d.runid=@rid and d.runid=s.runid and d.runid=w.runid
       order by d.node,d.loopnum,d.stage;
-- 62 rows affected (0.36 sec)
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'stage overhead' as name,d.id,
              d.duration-s.duration as duration
       from durations d, durations s
       where d.name='stage' and  s.name='app' and 
             d.loopnum=s.loopnum and d.stage=s.stage and d.node=s.node and
             (d.node!=10 or d.stage!=5) and d.stage!=1 and 
             d.runid=@rid and d.runid=s.runid
       order by d.node,d.loopnum,d.stage;
-- 1116 rows affected (2.05 sec)

-- mw overhead per loop:
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'loop overhead' as name,d.id,
              d.duration-s.duration-w.duration as duration
       from durations d, durations s, durations w
       where d.name='master loop' and  s.name='app loop' and 
             w.name='loop event wait' and 
             d.loopnum=s.loopnum and s.loopnum=w.loopnum and 
             d.node=s.node and s.node=w.node and 
             d.runid=@rid and d.runid=s.runid and d.runid=w.runid
       order by d.loopnum,d.stage,d.node;
-- 186 rows affected (2.51 sec)
