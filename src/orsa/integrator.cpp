#include <orsa/integrator.h>

#include <orsa/bodygroup.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

#include <iostream>

using namespace orsa;

Integrator::Integrator() : Referenced(true) {
    progressiveCleaningSteps = 0;
    onlyIfExtending = true;
}

bool Integrator::integrate(orsa::BodyGroup  * bg,
                           const orsa::Time & start,
                           const orsa::Time & stop,
                           const orsa::Time & sampling_period) {
    
    if (!bg->haveDynamicBodies(false)) {
        // no integration needed
        // ORSA_DEBUG("------ LEAVING EARLY -----");
        return true;
    }
    
    // #warning "remember: if any of the bodies has propulsion, care should be taken for timestep management,
    // to hit all the intervals where propulsion is used"
    
    doAbort = false;
    
    if (0) {
        // debug: all bodies in integration
        const orsa::BodyGroup::BodyList bl = bg->getBodyList();
        orsa::BodyGroup::BodyList::const_iterator it = bl.begin();
        while (it != bl.end()) {
            ORSA_DEBUG("integration includes body [%s] id: %i   time interval: [below]",
                       (*it)->getName().c_str(),
                       (*it)->id());
            orsa::print(start);
            orsa::print(stop);
            ++it;
        }
    }
    
    /* ORSA_DEBUG("CALL start: %f [day]",
       FromUnits(start.get_d(),Unit::DAY,-1));
       orsa::print(start);
       ORSA_DEBUG("CALL stop:  %f [day]",
       FromUnits(stop.get_d(),Unit::DAY,-1));
       orsa::print(stop);
    */
  
    // always start <= stop
    // but integration direction may not be from start to stop.... 
    // maybe I should change the name of the variables!!
    if (start > stop) {
        // ORSA_DEBUG("inverting...");
        return integrate(bg,stop,start,sampling_period);
    }
    
    // first, we need to have a common interval for the massive bodies
    if (bg->haveDynamicBodies(true)) {
        Time t1M, t2M;
        if (!(bg->getCommonInterval(t1M, t2M, true))) {
            ORSA_ERROR("cannot find common interval [massive bodies only]");
            return false;
        }
    }
    //
    // second, we need to have a common interval for all the bodies
    // if no interval is found, we must integrate massive bodies up to the epoch of the massless ones
    if (bg->haveDynamicBodies(false)) {
        bool needMicroIntegrations=false;
        {
            Time t1, t2;
            if (!(bg->getCommonInterval(t1, t2, false))) {
                needMicroIntegrations=true;
            }
        }
        if (needMicroIntegrations) {
            // compute global interval
            Time t1, t2;
            bg->getGlobalInterval(t1, t2, false);
            const orsa::Time t1G=t1;
            const orsa::Time t2G=t2;
            // extend them to include start,stop if necessary
            /* Time t1L=t1G;
               Time t2L=t2G;
               const Time t_min = std::min(start,stop);
               const Time t_max = std::max(start,stop);
               if (t_min < t1L) t1L=t_min;
               if (t_max > t2L) t2L=t_max;
            */
            
            // dump all massless bodies here
            osg::ref_ptr<orsa::BodyGroup> bgDump = new orsa::BodyGroup;
            {
                const orsa::BodyGroup::BodyList bl = bg->getBodyList();
                orsa::BodyGroup::BodyList::const_iterator it = bl.begin();
                while (it != bl.end()) {
                    if((*it)->getInitialConditions().inertial->mass()==0.0) {
                        bgDump->addBody((*it).get());
                        bg->removeBody((*it).get());
                    }
                    ++it;
                }
            }
            
            // integrate massive bodies alone (necessary?)
            /* if (!integrate(bg, t1***, t2***, sampling_period)) {
               ORSA_DEBUG("problems...");
               return false;
               }
            */
            
            // try to add massless bodies at small groups and integrate
            // the test for a good addition is that getCommonInterval returns true
            {
                // orsa::Cache<unsigned int> old_bgDump_size; // don't init now, to skip first iteration, where the massive bodies are integrated, along with possibly some massless bodies
                while (bgDump->size() != 0) {
                    const orsa::BodyGroup::BodyList bl = bgDump->getBodyList();
                    orsa::BodyGroup::BodyList::const_iterator it = bl.begin();
                    while (it != bl.end()) {
                        bg->addBody((*it).get());
                        bgDump->removeBody((*it).get());
                        Time t1, t2;
                        if (!(bg->getCommonInterval(t1, t2, false))) {
                            /* ORSA_DEBUG("cannot add body [%s] id: %i   [retrying later]",
                               (*it)->getName().c_str(),
                               (*it)->id());
                            */
                            // cannot add this body now, move it back into bgDump
                            bgDump->addBody((*it).get());
                            bg->removeBody((*it).get());
                        } else {
                            /* ORSA_DEBUG("added body [%s] id: %i",
                               (*it)->getName().c_str(),
                               (*it)->id());
                            */
                        }
                        ++it;
                    }
                    if (!integrate(bg, t1G, t2G, sampling_period)) {
                        ORSA_DEBUG("problems...");
                        return false;
                    }
                    /* if (old_bgDump_size.isSet()) {
                       if (bgDump->size() == old_bgDump_size.getRef()) {
                       ORSA_DEBUG("problem while dealing with massless bodies");
                       break;
                       }
                       }
                       old_bgDump_size = bgDump->size();
                    */
                }
            }
        }
    }
    
    // finally, we can retrieve the common interval for all bodies
    Time t1, t2;
    if (!(bg->getCommonInterval(t1, t2, false))) {
        ORSA_ERROR("cannot find common interval [all bodies]");
        return false;
    }
    //
    // orsa::print(t1);
    // orsa::print(t2);
    
    /* 
       if (t0 != start) {
       ORSA_DEBUG("splitting integration...");
       // two separate integrations: the first forward, the second backward in time
       // inverted order...
       return (integrate(bg,t0,start,sampling_period) && 
       integrate(bg,t0,stop, sampling_period));
       }
    */
    //
    bool backwardSplit = false;
    if ((start < t1) && (t1 != stop)) {
        // ORSA_DEBUG("splitting integration: start < t1");
        // integrate(bg,start,t1,sampling_period);
        backwardSplit = true;
    }
    //
    bool forwardSplit = false;
    if ((stop > t2) && (t2 != start)) {
        // ORSA_DEBUG("splitting integration: stop > t2");
        // integrate(bg,t2,stop, sampling_period);
        forwardSplit = true;
    }
    //
    if (backwardSplit) {
        if (forwardSplit) {
            // ORSA_DEBUG("B+F");
            return (integrate(bg,start,t1,sampling_period) &&
                    integrate(bg,t2,stop, sampling_period));
        } else {
            // ORSA_DEBUG("B");
            return integrate(bg,start,t1,sampling_period);
        }
    } else if (forwardSplit) {
        // ORSA_DEBUG("F");
        return integrate(bg,t2,stop, sampling_period);
    }
    
    Cache<mpz_class> sign;
    Time t;
    Time tEnd;
    //
    if (t1 == stop) {
        t    = t1;
        tEnd = start;
        sign.set(-1);
        // ORSA_DEBUG("--SIGN--MARK--");
    } else if (t2 == start) {
        t    = t2;
        tEnd = stop;
        sign.set(+1);
        // ORSA_DEBUG("--SIGN--MARK--");
    } else {
    
        // this is the case of a call for refinement, to create one more accurate snapshot of the system
        orsa::Time ctstart, ctstop;
        if (!bg->getClosestCommonTime(ctstart, start, false)) {
            ORSA_DEBUG("problems...");
        }
        if (!bg->getClosestCommonTime(ctstop, stop, false)) {
            ORSA_DEBUG("problems...");
        }
    
        /* 
           ORSA_DEBUG("ctstart: %f [day]",
           FromUnits(FromUnits(ctstart.getMuSec().get_d(),Unit::MICROSECOND),Unit::DAY,-1));
           ORSA_DEBUG(" ctstop: %f [day]",
           FromUnits(FromUnits(ctstop.getMuSec().get_d(),Unit::MICROSECOND),Unit::DAY,-1));
        */
    
        if ((ctstart==start) && (ctstop==stop)) {
            // ORSA_DEBUG("--- leaving ---");
            return true;
        } else if (ctstart==start) {
            // ORSA_DEBUG("new code ctstart...");
            t    = ctstop;
            tEnd = stop;
            sign.set((t-tEnd).getMuSec()/abs((t-tEnd).getMuSec()));
        } else if (ctstop==stop) {
            // ORSA_DEBUG("new code ctstop...");
            t    = ctstart;
            tEnd = start;
            sign.set((t-tEnd).getMuSec()/abs((t-tEnd).getMuSec()));
        } else {
            ORSA_DEBUG("---- CASE NOT HANDLED YET ----");
        }
    }	
  
    /* 
       mpz_class sign;
       //
       if (start < stop) {
       sign = 1;
       } else {
       sign = -1;
       }
    */
  
    const Time zeroTime = Time(0);
  
    // always start <= stop
    /* 
    // CANNOT DO THIS!!  
    if (start > stop) {
    ORSA_DEBUG("inverting...");
    return integrate(bg,stop,start,sampling_period);
    }
    */
  
    //
    /* 
       ORSA_DEBUG("t1: %f [day]",
       FromUnits(FromUnits(_t1.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1)());
       ORSA_DEBUG("t2: %f [day]",
       FromUnits(FromUnits(_t2.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1)());
    */
  
    /* 
       if ((start != _t1) || (_t1 < _t2)){
       ORSA_DEBUG("simple assumptions not working...");
       return false;
       }
    */
  
    // reset variables, needed when the same integrator is used multiple times
    // ORSA_DEBUG("calling Integrator::reset()...");
    // 
    reset();
  
    // integration statistics
    osg::ref_ptr< Statistic<double> > stat = new Statistic<double>;
    //
    osg::ref_ptr< RunningStatistic<double> > rs = new RunningStatistic<double>;
    rs->setLength(10);  
  
    // orsa::Time dt = sampling_period;
    //
    orsa::Time call_dt = sampling_period*sign.getRef();
    orsa::Time next_dt = sampling_period*sign.getRef();
  
    /* 
       if (start > stop) {
       call_dt = -call_dt;
       next_dt = -next_dt;
       }
    */
  
    bool ret_val = true;
  
    // const Time sampling_period(0,0,0,15,0);  // d,H,M,S,mu_S
 
    // while (t <= stop) {
  
    // a better first step
    if ((t+next_dt-tEnd)*sign.getRef() > zeroTime) {
        next_dt = tEnd - t;
    }
  
    while((tEnd-t)*sign.getRef() > zeroTime) {
    
        /* 
           ORSA_DEBUG("next_dt: %f   (mu: %Zi)",
           FromUnits(FromUnits(next_dt.getMuSec(),Unit::MICROSECOND),Unit::SECOND,-1)(),
           next_dt.getMuSec().get_mpz_t());
        */
    
        if (doAbort) {
            // ORSA_DEBUG("OUT//1");
            ret_val = false;
            break;
        }
    
        if (next_dt == zeroTime) {
            // ORSA_DEBUG("OUT//2");
            ret_val = true;
            break;
        }
    
        call_dt = next_dt;
    
        /* ORSA_DEBUG("call timestep: %f [day]     = %f",
           FromUnits(FromUnits(call_dt.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1)(),
           call_dt.get_d());
        */
    
        // while (_t <= (stop+sampling_period)) { // one more step...
        // std::cerr << "Integrator::integrate(...) --> t: " << FromUnits(FromUnits(_t.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1) << " [day]" << std::endl;
        
        // orsa::print(t);
        
        // if (!(step(bg,_t,dt))) {
        if (!(step(bg,t,call_dt,next_dt))) {
            ORSA_ERROR("problem encountered in step(...) call");
            // ORSA_DEBUG("OUT//3");
            ret_val = false;
            break;
        }
    
        if (0) {
            // debug
            BodyGroup::BodyList::const_iterator _b_it = bg->getBodyList().begin();
            while (_b_it != bg->getBodyList().end()) { 
                orsa::BodyGroup::BodyInterval * _b_interval = bg->getBodyInterval((*_b_it).get());
                orsa::BodyGroup::BodyInterval::DataType & _b_interval_data = _b_interval->getData();
                orsa::BodyGroup::BodyInterval::DataType::iterator _b_interval_data_it = _b_interval_data.begin();
                unsigned int counter=0;
                while (_b_interval_data_it != _b_interval_data.end()) {
                    ++counter;
                    ORSA_DEBUG("checking interval [%05i/%05i] [%s] time: %.6f",
                               counter,
                               (*_b_interval).size(),
                               (*_b_it)->getName().c_str(),
                               (*_b_interval_data_it).time.getRef().get_d());
                    ++_b_interval_data_it;
                }
                ++_b_it;	
            }
        }
    
        /* 
           ORSA_DEBUG("next timestep: %f [day]",
           FromUnits(FromUnits(next_dt.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1)());
        */
    
        // next_dt = call_dt;
    
        if (!lastCallRejected()) {
      
            singleStepDone(bg,t,call_dt,next_dt);
      
            t += call_dt;
      
            stat->insert(call_dt.get_d());
      
            // debug
            /* ORSA_DEBUG("t: %f [day]   E: %24.16Fe",
               orsa::FromUnits(t.get_d(),orsa::Unit::DAY,-1)(),
               bg->totalEnergy(t)());
            */

            // orsa::print(t);
            
            /* 
               ORSA_DEBUG("successful dt: %f +/- %f [day]",
               FromUnits(stat->average(),Unit::DAY,-1)(),
               FromUnits(stat->averageError(),Unit::DAY,-1)());
            */
            //
      
            if (progressiveCleaningSteps.isSet()) {
                if (progressiveCleaningSteps.getRef() > 0) {
                    if ((stat->entries() % progressiveCleaningSteps.getRef()) == 0) {
                        BodyGroup::BodyList::const_iterator _b_it = bg->getBodyList().begin();
                        while (_b_it != bg->getBodyList().end()) { 
                            if (!((*_b_it)->getInitialConditions().dynamic())) { 
                                ++_b_it;
                                continue;
                            }
                            orsa::BodyGroup::BodyInterval * _b_interval = bg->getBodyInterval((*_b_it).get());
                            orsa::BodyGroup::BodyInterval::DataType & _b_interval_data = _b_interval->getData();
                            orsa::BodyGroup::BodyInterval::DataType::iterator _b_interval_data_it = _b_interval_data.begin();
                            while (_b_interval_data_it != _b_interval_data.end()) {
                                /* ORSA_DEBUG("testing: body [%s]   t: %f   [tmp: %i]   this: %x   end: %x",
                                   (*_b_it)->getName().c_str(),
                                   (*_b_interval_data_it).time.getRef().get_d(),
                                   (*_b_interval_data_it).tmp,
                                   (&(*_b_interval_data_it)),
                                   (&(*(_b_interval_data.end()))));
                                */
                                if ((*_b_interval_data_it).tmp==true) {
                                    _b_interval_data_it = _b_interval_data.erase(_b_interval_data_it);
                                } else {
                                    ++_b_interval_data_it;
                                }
                            }     
                            // IMPORTANT!
                            _b_interval->update();
                            ++_b_it;
                        }
                    }
                }
            }

            if (keepOnlyLastStep.isSet()) {
                if (keepOnlyLastStep.getRef()) {
                    BodyGroup::BodyList::const_iterator _b_it = bg->getBodyList().begin();
                    while (_b_it != bg->getBodyList().end()) { 
                        if (!((*_b_it)->getInitialConditions().dynamic())) { 
                            ++_b_it;
                            continue;
                        }
                        orsa::BodyGroup::BodyInterval * _b_interval = bg->getBodyInterval((*_b_it).get());
                        orsa::BodyGroup::BodyInterval::DataType & _b_interval_data = _b_interval->getData();
                        orsa::BodyGroup::BodyInterval::DataType::iterator _b_interval_data_it = _b_interval_data.begin();
                        while (_b_interval_data_it != _b_interval_data.end()) {
                            /* ORSA_DEBUG("testing: body [%s]   t: %f   [tmp: %i]   this: %x   end: %x",
                               (*_b_it)->getName().c_str(),
                               (*_b_interval_data_it).time.getRef().get_d(),
                               (*_b_interval_data_it).tmp,
                               (&(*_b_interval_data_it)),
                               (&(*(_b_interval_data.end()))));
                            */
                            if ((*_b_interval_data_it).time.getRef()!=t) {
                                _b_interval_data_it = _b_interval_data.erase(_b_interval_data_it);
                                // ORSA_DEBUG("removed");
                            } else {
                                ++_b_interval_data_it;
                                // ORSA_DEBUG("kept...");
                            }
                        }     
                        // IMPORTANT!
                        _b_interval->update();
                        ++_b_it;
                    }	  
                }
            }
            
            rs->insert(call_dt.get_d());
            //
            /* if (rs->isFull()) {
               ORSA_DEBUG("radt: %f %f %f",
               FromUnits(t.get_d(),Unit::DAY,-1),
               FromUnits(rs->average(),Unit::SECOND,-1),
               FromUnits(rs->averageError(),Unit::SECOND,-1));
               }
            */
            
        } else {
            // last call was rejected...
        }
    
        if ((t+next_dt-tEnd)*sign.getRef() > zeroTime) {
            next_dt = tEnd - t;
            if (next_dt == zeroTime) {
                // out of here!
                // ORSA_DEBUG("OUT//4");
                break;
            }
        } else {
            orsa::Time eventTime;
            BodyGroup::BodyList::const_iterator bl_it = bg->getBodyList().begin();
            while (bl_it != bg->getBodyList().end()) { 
                if ((*bl_it)->propulsion.get()) {
                    eventTime = t;
                    if ((*bl_it)->propulsion->nextEventTime(eventTime,sign.getRef())) {
                        if ((t+next_dt-eventTime)*sign.getRef() > zeroTime) {
                            next_dt = eventTime - t;
                            //
                            // must call reset to allow step rejection
                            reset();
                        }
                    }
                }
                ++bl_it;
            }
        }
    }
    
    {
        BodyGroup::BodyList::const_iterator _b_it = bg->getBodyList().begin();
        while (_b_it != bg->getBodyList().end()) { 
            if (!((*_b_it)->getInitialConditions().dynamic())) { 
                ++_b_it;
                continue;
            }
            orsa::BodyGroup::BodyInterval * _b_interval = bg->getBodyInterval((*_b_it).get());
            orsa::BodyGroup::BodyInterval::DataType & _b_interval_data = _b_interval->getData();
            orsa::BodyGroup::BodyInterval::DataType::iterator _b_interval_data_it = _b_interval_data.begin();
            while (_b_interval_data_it != _b_interval_data.end()) {
                /* ORSA_DEBUG("testing: body [%s]   t: %f   [tmp: %i]   this: %x   end: %x",
                   (*_b_it)->getName().c_str(),
                   (*_b_interval_data_it).t.get_d(),
                   (*_b_interval_data_it).tmp,
                   (&(*_b_interval_data_it)),
                   (&(*(_b_interval_data.end()))));
                */
                if ((*_b_interval_data_it).tmp==true) {
                    _b_interval_data_it = _b_interval_data.erase(_b_interval_data_it);
                } else {
                    ++_b_interval_data_it;
                }
            }     
            // IMPORTANT!
            _b_interval->update();
            ++_b_it;
        }
    }
    
    // ORSA_DEBUG("--- done ---");
    return ret_val;
}
