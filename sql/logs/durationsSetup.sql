-- clean up from previous run
drop procedure fillloop;
drop procedure fillslices;
drop procedure fillallnodes;
drop procedure fixendloop;
drop view apploop;
drop view startloop;
drop view runtbl;


-- create the table
select min(micros) from logger where runid='rlp0130';
create table rlp0130 as select *, (micros-1201979215093257)/1000000.0 
       as secs from logger where runid='rlp0130';
alter table rlp0130 add column loopnum int(5), add column stage int(3);
create view runtbl as select * from rlp0130;

-- some indexes will help
create index nodeidx on rlp0130 (node);
create index sliceidx on rlp0130 (sliceid);

-- fill in the node column
update runtbl set 
       node=convert(substr(hostId from 5 for locate('.',hostId)-5), signed);

create view startloop as select * from runtbl 
                         where comment like '%ing stage loop number%';

delimiter //
create procedure fillloop(in nd int, in sl int)
begin
  declare cus varchar(1024);
  declare comm varchar(1024);
  declare ln int(5);
  declare stg int(3);
  declare loopst decimal(25,4);
  declare stagest decimal(25,4);
  declare t decimal(25,4);
  declare done int default 0;
  declare rec cursor for select secs,custom,comment from runtbl 
                         where node=nd and sliceid=sl
                         order by secs;
  declare continue handler for not found set done = 1;

  open rec;

  set ln=0;
  set stg = 0;
  set loopst = 0.0;
  set stagest = 0.0;
  repeat
    fetch rec into t,cus,comm;
    if cus like 'loop%' then
       update runtbl set loopnum=ln where node=nd and sliceid=sl and
                                          secs >= loopst and secs < t;
       set ln=convert(substr(cus from locate('||',cus)+2 
                                 for  locate('~~',cus)-(locate('||',cus)+2)),
                      signed);
       set loopst = t;
    end if;
    if comm like '%iStage %' then
       update runtbl set stage=stg where node=nd and sliceid=sl and
                                         secs >= stagest and secs < t;
       set stg=convert(substr(comm from locate('iStage ',comm)+7
                              for locate('~~',comm)-(locate('iStage ',comm)+7)),
                       signed);
       set stagest = t;
    end if;
  until done end repeat;
  update runtbl set loopnum=ln where node=nd and sliceid=sl and
                                      secs >= loopst and secs <= t;
  update runtbl set stage=stg where node=nd and sliceid=sl and
                                     secs >= stagest and secs <= t;
  close rec;
end;
//

create procedure fillslices(in nd int)
begin
  declare sl int(4);
  declare done int default 0;
  declare slrec cursor for select sliceid from runtbl where node=nd
                                  group by sliceid order by sliceid;
  declare continue handler for not found set done = 1;

  open slrec;

  repeat
    fetch slrec into sl;
    call fillloop(nd, sl);
  until done end repeat;

  close slrec;
end;
//

create procedure fillallnodes()
begin
  declare nd int(4);
  declare done int default 0;
  declare ndrec cursor for select node from runtbl where custom like 'loop%' 
                                  group by node order by node;
  declare continue handler for not found set done = 1;

  open ndrec;

  repeat 
    fetch ndrec into nd;
    call fillslices(nd);
  until done end repeat;

  close ndrec;
end;
//

create procedure fixendloop()
begin 
  declare eid int(11);
  declare lp int(5);
  declare done int default 0;
  declare endrec cursor for 
          select e.id,s.loopnum from startloop s, startloop e  
                 where s.sliceid=-1 and e.sliceid=-1 and s.node=e.node and
                       s.secs < e.secs and s.comment=e.comment;
  declare continue handler for not found set done = 1;

  open endrec;

  repeat
    fetch endrec into eid, lp;
    update runtbl set comment=concat('ending stage loop number ',
                                      format(lp,0), '~~')
           where id=eid;
  until done end repeat;

  close endrec;
end;
//

delimiter ;

CREATE TABLE `durations` (
  `runid` varchar(80) default NULL,
  `sliceid` int(11) default NULL,
  `date` varchar(30) default NULL,
  `micros` bigint(25) unsigned default NULL,
  `node` int(10) default NULL,
  `visitId` int(11) default NULL,
  `secs` decimal(25,4) default NULL,
  `loopnum` int(5) default NULL,
  `stage` int(3) default NULL,
  `name` varchar(64) default NULL,
  `id` int(11) default NULL,
  `duration` decimal(25,4) default NULL
);

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
       order by t1.node,t1.loopnum,t1.sliceid;

-- loop event wait:
create view loopeventwait as
       select node,loopnum,sum(duration) as sum from durations 
              where name='event wait' 
              group by node,loopnum order by node,loopnum;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'loop event wait' as name,d.id,s.sum as duration
       from durations d, loopeventwait s
       where d.stage=1 and d.loopnum=s.loopnum and d.node=s.node and 
             d.name='event wait'
       order by d.node,d.loopnum;

-- max process():
create view maxproc as 
       select loopnum,stage,max(duration) as maxdur from durations 
       where name='stage process' and node!=9 and node!=10 
       group by loopnum,stage;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'max stage process' as name,d.id,d.duration 
       from durations d, maxproc m 
       where d.loopnum=m.loopnum and d.stage=m.stage and 
             d.duration=m.maxdur and d.name='stage process' and 
             d.micros!=1201849713848155 and d.micros!=1201851085989183 and
             d.micros!=1201851085988935
       order by loopnum,stage;
create view maxproc9 as 
       select loopnum,stage,max(duration) as maxdur from durations 
       where name='stage process' and node=9
       group by loopnum,stage;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'max stage process' as name,d.id,d.duration 
       from durations d, maxproc9 m 
       where d.loopnum=m.loopnum and d.stage=m.stage and d.node=9 and 
             d.duration=m.maxdur and d.name='stage process' 
       order by loopnum,stage;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'max stage process' as name,d.id,d.duration 
       from durations d 
       where d.name='stage process' and d.node=10
       order by loopnum,stage;

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
             d.node=s.node and s.node!=9 and s.node!=10 and 
             p.node!=9 and p.node!=10 
       order by d.loopnum,d.stage,d.node;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'app' as name,d.id,
              d.duration+s.duration+p.duration as duration
       from durations d, durations s, durations p 
       where d.name='stage preprocess' and  s.name='stage postprocess' and 
             p.name='max stage process' and 
             d.loopnum=s.loopnum and s.loopnum=p.loopnum and 
             d.stage=s.stage and s.stage=p.stage and 
             d.node=s.node and s.node=p.node and s.node=9
       order by d.loopnum,d.stage,d.node;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'app' as name,d.id,
              d.duration+s.duration+p.duration as duration
       from durations d, durations s, durations p 
       where d.name='stage preprocess' and  s.name='stage postprocess' and 
             p.name='max stage process' and 
             d.loopnum=s.loopnum and s.loopnum=p.loopnum and 
             d.stage=s.stage and s.stage=p.stage and 
             d.node=s.node and s.node=p.node and s.node=10
       order by d.loopnum,d.stage,d.node;

--- app code over entire loop:
create view apploop as 
       select runid,node,loopnum,sum(duration) as sum from durations 
       where name='app' group by runid,node,loopnum 
       order by node,loopnum;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'app loop' as name,d.id,s.sum as duration
       from durations d, apploop s
       where d.stage=1 and d.loopnum=s.loopnum and d.node=s.node and 
             d.name='app'
       order by d.node,d.loopnum;

-- mw overhead per stage:
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'stage overhead' as name,d.id,
              d.duration-s.duration-w.duration as duration
       from durations d, durations s, durations w
       where d.name='stage' and  s.name='app' and w.name='event wait' and 
             d.loopnum=s.loopnum and s.loopnum=w.loopnum and 
             d.stage=1 and d.stage=s.stage and s.stage=w.stage and 
             d.node=s.node and s.node=w.node 
       order by d.node,d.loopnum,d.stage;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'stage overhead' as name,d.id,
              d.duration-s.duration-w.duration as duration
       from durations d, durations s, durations w
       where d.name='stage' and  s.name='app' and w.name='event wait' and 
             d.loopnum=s.loopnum and s.loopnum=w.loopnum and 
             d.stage=5 and d.stage=s.stage and s.stage=w.stage and 
             d.node=10 and d.node=s.node and s.node=w.node 
       order by d.node,d.loopnum,d.stage;
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'stage overhead' as name,d.id,
              d.duration-s.duration as duration
       from durations d, durations s
       where d.name='stage' and  s.name='app' and 
             d.loopnum=s.loopnum and d.stage=s.stage and d.node=s.node and
             (d.node!=10 or d.stage!=5) and d.stage!=1 
       order by d.node,d.loopnum,d.stage;

-- mw overhead per loop:
insert into durations
       select d.runid,d.sliceid,d.date,d.micros,d.node,d.visitId,d.secs,
              d.loopnum,d.stage,'loop overhead' as name,d.id,
              d.duration-s.duration-w.duration as duration
       from durations d, durations s, durations w
       where d.name='master loop' and  s.name='app loop' and 
             w.name='loop event wait' and 
             d.loopnum=s.loopnum and s.loopnum=w.loopnum and 
             d.node=s.node and s.node=w.node
       order by d.loopnum,d.stage,d.node;
