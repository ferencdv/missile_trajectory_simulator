"""The numerical simulation. Basic text interface provided when run as main. Real interface in gui.pyw"""

from math import *
global VBO
global calcrange
global burnout_angle
global opt_burnout_angle
global m_old

class Simulation(object):
    """The numerical simulation"""
    def __setattr__(self,name,value):
        object.__setattr__(self, name, value)
        #used to allow attribute write access by other methods
    
    def __init__(self,parent):
        self.parent = parent #for reference to gui appframe
        #stage parameter lists
        self.burntime = ['']
        self.thrust0 = ['']
        self.Isp0 = ['']
        self.m0 = ['']
        self.fuelfraction = ['']
        self.fuelmass = ['']
        self.dMdt = ['']
        #fins parameter lists
        self.FinHeight = ['']
        self.SweepAngle = ['']
        self.FinTipChordLength = ['']
        self.FinRootChord = ['']
        self.FinThickness = ['']
        #results dict
        self.data = {'Time':[0],'Height':[0],'Mass':[0],'Velocity':[0],'Thrust':[0],'Drag':[0],'CD':[0],'Gamma':[pi/2],'Range':[0]}
        self.Nosecone = ['']
    def integrate(self,trajectory):
        t = 0.0     # time
        v = 0       # initial v
        h = 0.001   # initial h must be small but non-zero
        psi = 0     # range angle: range = psi * Rearth
        rho = 0.0   # air density at current altitude
        p_height = 0.0 # air pressure at current altitude
        gamma = self.to_radians(90) #launch angle, from horizontal
        
        
        #print "Start Simulation"
        # ref for printing results in GUI mode
        try:
            app = wx.GetTopLevelParent(self.parent)
        except NameError:
            pass
        self.trajectory = trajectory #make ref for eta function
        ##### SET INTEGRATION PARAMETERS
        tEND = 20000        #timeout value
        dtprint = 1         #time interval between printing output
        Htrans = 20000  #height [m] at which transition from laminar to turbulent heating occurs
        deltaend = .1       #time increment used for integration
        deltatinit = .01    #time increment for t < tinit + 1 sec
        mtot = 0.0
        burntimetot = 0.0
        tinit = burntimetot + 1 # integrate more carefully during burn
        #####
        apogee = 0.0
        Thrust = 0.0
        drag = 0.0
        cd = 0.0
        ##### SET CONSTANTS
        Rearth = 6370000 #[m]
        g0 = 9.8066 #[m/s^2]
        #
        ##### INITIALIZE ROCKET MODEL
                
        for i in range(1,self.numstages+1):
            print 'Inside initialize rocket model=',self.numstages,i
            print '>=>',i,self.m0[i],self.dMdt[i],self.thrust0[i],self.fuelmass[i]
            mtot += self.m0[i] #sum total mass
            self.burntime.append(self.Isp0[i]*9.81*self.fuelmass[i]/self.thrust0[i])
            burntimetot += self.burntime[i] #sum total burn time
        mtot += self.payload

        ##### FIN DATA
        fh = self.FinHeight # this is S(m)
        sa = self.SweepAngle # Angle in degrees
        CTm = self.FinTipChordLength # CT[m]
        CRm = self.FinRootChord # CR[m]
        TFm = self.FinThickness # TF[m]
        #numfin = self.FinsChoice # number of fins
        
        global fh, sa, CTm, CRm, TFm
        print "FH=",fh,sa,CTm,CRm,TFm


        area_missile = (self.missilediam/2)**2 * pi #[m^2]
        area_rv = self.rvdiam/2**2 * pi #[m^2]
        ##### CALCULATE THE BURNOUT VELOCITY
        print 'Num Stages=',self.numstages,mtot

        if self.numstages == 1:
           VBO = self.Isp0[1]*9.81*log(mtot/(mtot - self.fuelmass[1]))
           print VBO
        if self.numstages ==2:
           VBO = self.Isp0[1]*9.81*log(mtot/(mtot - self.fuelmass[1]))
           VBO = VBO + self.Isp0[2]*9.81*log((mtot-self.m0[1])/(mtot - self.m0[1] - self.fuelmass[2]))
           print VBO
        if self.numstages ==3:
           VBO = 0.
           VBO = self.Isp0[1]*9.81*log(mtot/(mtot - self.fuelmass[1]))
           VBO = VBO + self.Isp0[2]*9.81*log((mtot-self.m0[1])/(mtot - self.m0[1] - self.fuelmass[2]))
           VBO = VBO + self.Isp0[3]*9.81*log((mtot-self.m0[1]-self.m0[2])/(mtot -self.m0[1] -self.m0[2] - self.fuelmass[3]))
           print VBO
        if self.numstages ==4:
           VBO = 0.
           VBO = self.Isp0[1]*9.81*log(mtot/(mtot - self.fuelmass[1]))
           VBO = VBO + self.Isp0[2]*9.81*log((mtot-self.m0[1])/(mtot - self.m0[1] - self.fuelmass[2]))
           VBO = VBO + self.Isp0[3]*9.81*log((mtot-self.m0[1]-self.m0[2])/(mtot -self.m0[1] -self.m0[2] - self.fuelmass[3]))
           VBO = VBO + self.Isp0[4]*9.81*log((mtot-self.m0[1]-self.m0[2]-self.m0[3])/(mtot -self.m0[1] -self.m0[2] - self.m0[3]-self.fuelmass[4]))
           print VBO
        if self.numstages ==5:
           VBO = 0.
           VBO = self.Isp0[1]*9.81*log(mtot/(mtot - self.fuelmass[1]))
           VBO = VBO + self.Isp0[2]*9.81*log((mtot-self.m0[1])/(mtot - self.m0[1] - self.fuelmass[2]))
           VBO = VBO + self.Isp0[3]*9.81*log((mtot-self.m0[1]-self.m0[2])/(mtot -self.m0[1] -self.m0[2] - self.fuelmass[3]))
           VBO = VBO + self.Isp0[4]*9.81*log((mtot-self.m0[1]-self.m0[2]-self.m0[3])/(mtot -self.m0[1] -self.m0[2] - self.m0[3]-self.fuelmass[4]))
           VBO = VBO + self.Isp0[5]*9.81*log((mtot-self.m0[1]-self.m0[2]-self.m0[3]-self.m0[4])/(mtot -self.m0[1] -self.m0[2] - self.m0[3]-self.m0[4]-self.fuelmass[5]))
           print 'VBO=',VBO/1000.
        # Now no need to estimate the range anymore and can remove est range button
        ##### INTEGRATE
        #   
        #Initialize variables
        deltat = deltatinit
        flagdeltat = True
        m = mtot
        print 'm=mtot',m

        #
        print 'dmdt=',self.dMdt[1]
        dMdt0 = self.dMdt[1]
        tprint = dtprint #tprint is time at which printing of output will next occur
        flag = True # controls printing parameters at burnout of stages
        tlimit = self.burntime[1] # ditto
        nstage = 1  # used at burnout of stages
        gamma_half = gamma # angle of missile or RV w/ local horizon
        
        if self.trajectory == 'Minimum Energy':
           #set burnout angle to optimum for MET
           #uses Wheelon's form of the equations
           calcrange = exp((2500.*(VBO/1000.)+23629.)/4477.)
           print 'calcrange=', calcrange 
           opt_burnout_angle = pi/2 - .25*(1000*calcrange/Rearth + pi)
           print 'opt burnout angle in MET',opt_burnout_angle*180/3.141592

        if self.trajectory == 'Burnout Angle':
           opt_burnout_angle = self.burnout_angle*3.14592/180.
           print 'opt burnout angle in BOA',opt_burnout_angle*180/3.141592
           #print 'burn out angle', burnout_angle
           #use this optimum burnout angle to linearize turn angle, from horizontal

        
        #Integrate
        while t < tEND and h > 0: # big loop
            #save data to Results dict
            self.data['Time'].append(t) #in tenths seconds
            self.data['Height'].append(h) #in meters
            self.data['Mass'].append(m) #in kg
            self.data['Velocity'].append(v) #in meters/second
            self.data['Thrust'].append(Thrust) #in in kgf
            self.data['Drag'].append(drag) #in N
            self.data['CD'].append(cd) # dimensionless
            #self.data['cdfins'].append(cdfins) # dimensionless
            self.data['Gamma'].append(gamma) #in degrees from horizontal
            self.data['Range'].append(Rearth*psi) #in meters
            
            if (t + deltat/5) >= tinit and flagdeltat == True:
                deltat = deltaend
                flagdeltat = False
             #
            # save old values
            psi_old = psi
            h_old = h
            gamma_old = gamma
            v_old = v
            m_old = m
            t_old = t
            #
            if (t + deltat/5) <= burntimetot: 
                m_half = m_old - (dMdt0 * deltat/2) #burn fuel
                area = area_missile
            else:
                area = area_rv
            #calculate drag and include cd
            rho = self.density(h)
			
			
			
            cd = self.Cdrag(v_old,h)
            drag = cd*area*rho*(v_old**2)/2
            
            """ # calculate thrust as function of altitude
            #NEW EQUATIONS, from Charles Vick
            h_vacuum = 160934 #~100 miles
            Thrust_ideal = self.Isp0[nstage]*self.dMdt[nstage]*9.81
            if (t + deltat/5) > burntimetot:
                Thrust_pct_increase = 0
                #out of fuel, no thrust
            elif h < h_vacuum:
                h_norm = h / h_vacuum
                Thrust_pct_increase = -.4339*(h_norm)**3+.6233*(h_norm)**2-.01*(h_norm)+1.004
                #3rd order polynomial line fit from Saturn-V data on thrust vs. height
                
            elif h > h_vacuum and nstage == 1:
                Thrust_pct_increase = 1.19
                Thrust = Thrust_ideal*Thrust_pct_increase
            elif nstage > 1:
                Thrust_pct_increase = 1
                #assuming that stage Isp is correct for vacuum
            #Thrust = Thrust_ideal*Thrust_pct_increase
            #Force = Thrust - drag
            #note that Force will be negative during reentry """
                     
			# Because this routine comes second I assume it reclaculates the thrust 
            #OLD EQUATIONS, from David Wright
            #requires us to know nozzle area, which we don't $ FDV: We actually do often so that is why I have included this but leave in the other code.
            p0 = self.pressure(0)
            p_height = self.pressure(h)
            nozarea=self.nozzlearea
            if (t + deltat/5) > burntimetot:
               Thrust = 0.0
            elif nstage == 1:  
               Thrust = self.Isp0[1]*self.dMdt[1]*9.81 + nozarea*pi*(p0-p_height)
            elif nstage > 1:
               Thrust = self.Isp0[nstage]*self.dMdt[nstage]*9.81
                                       
            Force = Thrust - drag
            #print 'STAGE=>',t,nstage,'ISP,dMdT',self.Isp0[nstage],self.dMdt[nstage],'Thrust:',Thrust,'Force:',Force,'v_old',v_old
            #note that Force will be negative during reentry
                
            #
            g = g0*Rearth**2/(h+Rearth)**2 #calculate grav accel at height
            
            ETA_old = self.eta(h_old,t_old)
            #
            # Integration is variant of Runge-Kutta-2.
            # 1- Calculate values at midpoint, t = t_old + deltat/2
            #
            t_half = t_old + deltat/2
            d_psi = (v_old * cos(gamma_old)/(Rearth + h_old)) * deltat/2
            psi_half = psi_old + d_psi
            h_half = h_old + v_old*sin(gamma_old)*deltat/2
            #
            # calculate gamma
            
            vertical_flight_period = 5 # seconds this is very arbitrary! REVISIT THIS
            if t < vertical_flight_period:
                #force gamma to be constant early in flight
                dgamma = 0.0
            elif (t >= vertical_flight_period) and (t <= burntimetot):
                dgamma = ((opt_burnout_angle - pi/2)/(burntimetot - vertical_flight_period))
                # Flying the missile at the burnout angle until end of burnout
            else:
                dgamma = d_psi/(deltat/2) + Force*sin(ETA_old)/(v_old * m_old) - (g*cos(gamma_old)/v_old)
            
            #integrate it
            gamma_half = gamma_old + dgamma*deltat/2
            
            # calculate dv
            dv = (Force/m_old)*cos(ETA_old) - g*sin(gamma_old)
                        
            v_half = v_old + dv*deltat/2
            #
            #
            # 2- Use derivatives at midpoint to calculate values at t + deltat
            ETA_half = self.eta(h_half,t_half)
            # Increment time
            t += deltat
            #
            d_psi_half = (v_half*cos(gamma_half))/(Rearth+h_half) * deltat
            psi = psi_old + d_psi_half
            h = h_old + v_half*sin(gamma_half)*deltat
            if h > h_old:
                apogee = h
                v_apogee = v

            vertical_flight_period = 5
            if t <= vertical_flight_period:
                dgamma_half = 0.0
            elif (t > vertical_flight_period) and (t <= burntimetot):
                global opt_burnout_angle
                dgamma_half = ((opt_burnout_angle - pi/2)/(burntimetot - vertical_flight_period))
            else:
                #use Wright's equation, hopefully not too disjoint with previous
                dgamma_half = d_psi_half/(deltat) + (Force/(v_half*m_half))*sin(ETA_half) - (g*cos(gamma_half)/v_half)
                
            gamma = gamma_old + dgamma_half*deltat

            if (t + deltat/5) <= burntimetot:
                m = m_old - dMdt0 * deltat
                #burn fuel mass 
    
            dv_half = (Force/m_half)*cos(ETA_half) - g*sin(gamma_half)
            v = v_old + dv_half*deltat
                        
            #Print data at stage burnout
            if (t + deltat / 5) > tlimit and flag == True:
                if __name__ == "__main__":
                    #Simple text printout
                    print "Stage %i burnout" % nstage
                    print "Velocity (km/s): ",v/1000
                    print "Angle (deg h): ",gamma*180/pi
                    print "Range (km): ",Rearth*psi/1000
                    print "Time (sec): ",t
                else:
                    #GUI printout
                    app.Results.StageVelocityResult[nstage].SetValue("%4.2f" % float(v/1000))
                    app.Results.StageAngleResult[nstage].SetValue("%4.2f" % float(gamma*180/pi))
                    app.Results.StageHeightResult[nstage].SetValue("%4.2f" % float(h/1000))
                    app.Results.StageRangeResult[nstage].SetValue("%4.2f" % float(Rearth*psi/1000))
                    app.Results.StageTimeResult[nstage].SetValue("%4.2f" % t)
 
                # m = mtot - self.m0[nstage]
                if nstage == 1:
                    m = mtot - self.m0[1]
                    print 'NSTAGE=',nstage, m
                if nstage == 2:
                    m = mtot - self.m0[1] - self.m0[2]
                    print 'NSTAGE=',nstage, m                    
                if nstage == 3:
                    m = mtot - self.m0[1] - self.m0[2] - self.m0[3]
                    print 'NSTAGE=',nstage, m                    
                m_old = m


                if nstage < self.numstages:
                    nstage += 1
                    tlimit += self.burntime[nstage] #set time to next print burnout
                    dMdt0 = self.dMdt[nstage]
                else:
                    flag = False
                
            #END BIG LOOP




        if t >= tEND:
            if __name__ == "__main__":
                print "Simulation exceeded time limit."
            else:
                dlg = wx.MessageDialog(self.parent,"Exceeded time limit, results are likely invalid.","Simulation error",wx.OK | wx.ICON_INFORMATION)
                dlg.ShowModal()
                dlg.Destroy()


        #print "Done"
        if __name__ == "__main__":
            #print final results
            print "Range (km): ",psi*Rearth/1000
            print "Apogee (km): ",apogee/1000
            print "Time to target (sec): ",t
        else:
            #put results in frame
            app.Results.ApogeeResult.SetValue("%4.2f" % float(apogee/1000))
            app.Results.ApogeeVelocityResult.SetValue("%4.3f" % float(v/1000))
            app.Results.RangeResult.SetValue("%4.3f" % float(Rearth*psi/1000))
            app.Results.FlightTimeResult.SetValue("%4.1f" % t)
            
        return (self.data)
            
                
    def eta(self,h,t):
         #for 11,000km MET, from <Gronlund and Wright, "Depressed Trajectory SLBMS", Science and Global Security, 1992, Vol 3, p101-159>
        # only used for Thrust Vector trajectories
        if self.trajectory == 'Thrust Vector':
            if t > self.TStartTurn and t < self.TEndTurn:
                eta = -self.to_radians(self.TurnAngle)
            else:
                eta = 0.0
        else:
            eta = 0.0
        return eta      
    
    def density(self,h):
        "Calculates air density at altitude"    
        rho0 = 1.225 #[kg/m^3] air density at sea level
        if h < 19200:
            #use barometric formula, where 8420 is effective height of atmosphere [m]
            rho = rho0 * exp(-h/8420)
        elif h > 19200 and h < 47000:
            #use 1976 Standard Atmosphere model
            #http://modelweb.gsfc.nasa.gov/atmos/us_standard.html
            #from http://scipp.ucsc.edu/outreach/balloon/glost/environment3.html
            rho = rho0 * (.857003 + h/57947)**-13.201
        else:
            #vacuum
            rho = 1.e-6
        return rho
        
    def temperature(self,h):
        "Calculates air temperature [Celsius] at altitude [m]"
        #from equations at 
        #   http://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
        if h <= 11000:
            #troposphere
            t = 15.04 - .00649*h
        elif h <= 25000:
            #lower stratosphere
            t = -56.46
        elif h > 25000:
            t = -131.21 + .00299*h
        return t
    
    def pressure(self,h):
        "Calculates air pressure [Pa] at altitude [m]"
        #from equations at 
        #   http://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
        
        t = self.temperature(h)
        
        if h <= 11000:
            #troposphere
            p = 101.29 * ((t+273.1)/288.08)**5.256
        elif h <= 25000:
            #lower stratosphere
            p = 22.65*exp(1.73-.000157*h)
        elif h > 25000:
            p = 2.488 * ((t+273.1)/288.08)**-11.388
        return p
        
        
        """ 
        Changes done on 5/13 base code was changed from gui-7-new-4.23.py 
        Note that this is a fit of total drag which includes the body drag
        plus the fricton drag and base drag. It does not include the fin drag.   
        Note that this drag also does not include the drag due to various length
        bodies because it does not seem to be a large effect.
         
        FITTING THE DATA
        
        The drag functions for various L/D (length of nosecone/ max diameter of nosecone) 
        and for various mach speeds were produced using HyperCFD code. These curves can
        be divided into three segments from M=0 to M=1.05 [1], from M=1.05 to M=1.2 [2]
        and from M=1.2 and above. In each of these sections the curves vary as a function
        of M and of L/D ratio. Segment [1] is the easiest and is interpolated linearly
        from L/D=1 to L/D=3. The next segment [2] is a very steep fall off to M=1.2 and 
        is also intepreted linearly from M=1.05 to M=1.2 (essentially a linear interpolation
        beteen these two points. 
        M=0  1.05   1.2                                                                                 20
        | 1  |  2   |   3                                                                                |
        ................................................................................................   
        .y````````:.````````+````````.+````````+````````./````````+`````````+````````:-````````+`````````/  
        .y        :.        +        `/        /        `/        /         /        :.        /         /  
        .y        :.        +        `/        /        `/        /         /        :.        /         /  
        .y        :.        +        `/        /        `/        /         /        :.        /         /  
        .h......../:........o........-o........o........-+........o.........o........+:........o........./  
        .y````````:.````````+`````````+````````+````````./````````+`````````+````````:-````````+`````````/  
        -y        :.        +        `/        /        `/        /         /        :.        /         /  
        -y .-.    :.        +        `/        /        `/        /         /        :.        /         /  
        -y /h/    :.        +        `/        /        `/        /         /        :.        /         /  
        -h--h----.+:........o........-o........o........:+........o.........o........+:.......-o........-/  
        -y``y`````:.````````+`````````+````````+````````./````````+`````````+````````:-````````+`````````/  
        -y  y`    :.        +        `/        /        `/        /         /        :.        /         /  
        -y  y`    :.        +        `/        /        `/        /         /        :.        /         /  
        -y  y`    :.        +        `/        /        `/        /         /        :.        /         /  
        -h--h/----+:--------s---------o--------o--------:o--------o---------o--------+/--------o---------/  
        -y  h:    :.        +        `/        /        `/        /         /        :.        /         /  
        -y  h/    :.        +        `/        /        `/        /         /        :.        /         /  
        -y  hy    :.        +        `/        /        `/        /         /        :.        /         /  
        -y `hy    :.        +        `/        /        ./        /         /        :.        /         /  
        -h--hh----+:--------o---------o--------o--------:+--------o---------o--------+:--------o---------/  
        -y  ys.   :.        +        `/        /        `/        /         /        :.        /         /  
        -y  ys-   :.        +        `/        /        `/        /         /        :.        /         /  
        -y  y-o   :.        +        `/        /        `/        /         /        :.        /         /  
        .y``y`o-``:.`      `+`      ``/      ` +``      ./       `+        `+`       :-       `/        `/  
        -h--h-+o--+:--------o---------o--------o--------:o--------o---------o--------+:--------o---------/  
        -y  y``s. :.        +        `/        /        `/        /         /        :.        /         /  
        -y  y` -o :.        +        `/        /        `/        /         /        :.        /         /  
        -y  y` .s`:.        +        `/        /        `/        /         /        :.        /         /  
        -m+.y...+//-........o.........+........+........-+........+.........+......../:........+........./  
        .m+`y.```os-````````o````````.+````````+````````./````````+`````````+````````/-````````+`````````/  
        .Nssy-   .d.        +        `/        /        `/        /         /        :.        /         /  
        -h.--`    ++        +        `/        /        `/        /         /        :.        /         /  
        -y        :h.       +        `/        /        `/        /         /        :.        /         /  
        .h--------++s-------o---------o--------o--------:o--------o---------o--------+:--------o---------/  
        .y````````:.-o.`````+`````````+````````+````````./````````+`````````+````````:-````````/`````````/  
        .y        :. -o`    +        `/        /        `/        /         /        :.        /         /  
        .y        :.  .o-   +        `/        /        `/        /         /        :.        /         /  
        .y        :.   `/+. +        `/        /        `/        /         /        :.        /         /  
        .h......../:....-:o+s-.......-o........o........:+........o.........o........+:........o........-/  
        .y````````:-```````-y+/-`````.+````````+````````./````````+`````````+````````/-````````+`````````/  
        .y        :.        +`-//+--../        /        `/        /         /        :.        /         /  
        .y        :.        +    `--:/s/:--..``+        `/        /         /        :.        /         /  
        .y        :.        +        `+`..---//s///::::::+---.....+..```````+```     :.        /         /  
        .h--------+:--------s---------o--------o---::::/+s///+++++y++///////s///+////o+////////s/////++++/  
        .y        :.        +        `/        /        `/        /         /        :-       `/    `..../  
        .y        :.        +        `/        /        `/        /         /        :.        /         /  
        .y        :.        +        `/        /        `/        /         /        :.        /         /  
        .h......../-........o........-+........+........-+........+.........+......../:........+........./  
        `/::::::::::::::::::/:::::::::/::::::::/::::::::::::::::::/:::::::::/::::::::::::::::::::::::::::. 

        
        These curves tend to have a very steep fall off from M=1 to 
        M=1.2 so it is not possible to have a good fit from M>1.05(ie [2] and [3]). The final
        fit is in segment [3] from M=1.2 onwards. The fit that is done generally is the Hoerl
        function which is:Cd = a*b^x*x^c, where a, b, c are all constants that vary with L/D.
        So the fit is done for each L/D = 1,1.5,2,2.5,3 for M=1.05,1.2,1.3,1.4,1.5,2,3,4,8,15 
        to a Hoerl function. Then the variation in the coefficients a, b, and c are fit which 
        relate the variation of the coefficients in L/D to the variation over speed. The fits
        that work vary for different nosecone models but it tends to fit:
        For a the function a*exp(x^b), for b simply an average (a constant), for c the logistic 
        formula which is 1/(1+b*exp(-c*x)), where in all cases x = L/D. Similarly, the linear 
        fits to M=1.05 to M=1.2 also give two coefficient m and b which are also fit to 
        reciprical functions ie. a = 1/(a+b*(L/D). For section [1] the fit between 0. and 1.05 
        is assumed to be a constant (generally between 0.05 and 0.2) that varies slightly between
        L/D=1 to L/D=3. We assume this to be linear so that the constant varies according to a
        linear function of CD at L/D at 1 and 3.

        The fits are done for 5 functions: Conic, Parabola, Elliptical, Sears-Haack and Tangent Ogive.

        What is not included yet is the drag due to the fins, which is a sizable effect and will
        be done next. 

        """

    def Cdrag (self,v,h):

        if v == 0.:
            v=0.0001
        t = self.temperature(h) + 273.15 #convert to kelvin
        a = sqrt(1.4*287*t) # 1.4 is the ratio of specific heats. See: https://www.grc.nasa.gov/WWW/BGH/specheat.html
        # 287 is the gas constant in air in J/kg/kelvin. See: https://www.grc.nasa.gov/WWW/BGH/eqstat.html
        mach = v/a
        ld=self.LdivD
        diam = self.missilediam       
        Rho = self.density(h)  
        #numfin = self.numfins
        
        #Drag function for V2
        #derived from Sutton, "Rocket Propulsion Elements", 7th ed, p108
        #probably not that relevant to other body types

        cdfins = 0.
        #print 'fh=',fh
        lflm = (fh**2+(0.5*(CTm-CRm)+fh/tan(3.141592*sa/180))**2)**(0.5)
        #print 'lfm=',lflm
		
        AFEP = 0.5*(CRm+CTm)*fh
        #print 'AFEP=',AFEP,CRm,CTm,sa
        AFP = AFEP + 0.5*CRm*diam
        Mu = 0.00001789 # kg/m.s
        #print 'Rho,a',Rho,v,lflm,Mu,a
        Renolds = Rho*v*lflm/Mu # Calculation for Renolds Numbr
        #print 'Renolds', Renolds
        RenoldsCrit= 500000 # Five Hundred Thousand
        FinBt= Renolds**(1./5.)
        FinB = RenoldsCrit*((0.074/(Renolds**(1./5.)))-(1.328/sqrt(Renolds)))
        # print 'Renoldscrit=',Renolds,(0.074/(FinBt)),((1.328/sqrt(Renolds))),FinB,FinBt
        CFforRenoldsLow = 1.328/sqrt(Renolds)
        CFforRenoldsHigh = (0.074/(FinBt))-FinB/Renolds
 
        if Renolds < RenoldsCrit:
            CF = CFforRenoldsLow
        else:
            CF = CFforRenoldsHigh   
        #print 'AFE=',AFEP,AFP,diam     
        cdfins = 2.*CF*(1.+2.*(TFm/lflm))*(4.*1.*(AFP))/(3.141592*diam**2.)

        #print 'AFE=',AFEP
        #print 'Renolds', Renolds
        #print 'CFforRenoldsLow',CFforRenoldsLow,CFforRenoldsHigh, Renolds
        #print 'Rho,Renolds',Rho,Renolds,'     FinB,CFforRenoldsLow,CFforRenoldsHigh,CF,AFP',FinB,CFforRenoldsLow,CFforRenoldsHigh,CF,AFP

        cdfins= 0.
      
        if self.Nosecone == 'V2':
            if mach > 5:
               cd = 0.15-cdfins
            elif mach > 1.8 and mach <= 5:
               cd = -0.03125*mach + 0.30625-cdfins
            elif mach > 1.2 and mach <= 1.8:
               cd = -0.25*mach + 0.7-cdfins
            elif mach > 0.8 and mach <= 1.2:
               cd = 0.625*mach - 0.35-cdfins
            elif mach <= 0.8:
               cd = 0.15-cdfins
			   
        if self.Nosecone == 'elliptical':    
            if mach >= 1.2: # All values taken from fits to HyperCFD
               CALDL = 0.824584774*ld**-0.532619017
               CBLDL = 1.0156845
               CCLDL = -0.226354-0.238389*log(ld)
               cd = (CALDL*CBLDL**mach)*mach**(CCLDL)
            elif mach >= 1.05 and mach <= 1.2:
               MCONIC = 1/(-0.2383263138-0.266070229318*ld)
               BCONIC = 1/(0.15266+0.160535*ld)
               cd = MCONIC*mach+BCONIC
            elif mach <1.05:
               cd = -0.05*ld + 0.25

        if self.Nosecone == 'Conical':   
            if mach > 1.5: # All values taken from fits to HyperCFD   
               ACONIC = 1.619038033*exp(-1.31926217*ld)
               BCONIC = ld/(-0.45318-0.89392*ld)
               CCONIC = 0.886118*exp(-ld/1.121185)
               cd = ACONIC/(1+BCONIC*exp(-CCONIC*mach))
            elif mach >= 1.05 and mach <= 1.5:
               MCONIC = 1/(-0.10823-0.81349*ld)
               BCONIC = 1/(0.054882+0.363845*ld)
               cd = MCONIC*mach+BCONIC
            elif mach <1.05:
               cd = 0.075*ld + 0.275

        if (self.Nosecone == 'tangent ogive'): # assume that both nosecone shapes give same shape
            if mach >= 1.2: # All values taken from fits to HyperCFD
               APARA = 0.278184983*exp(ld**-0.8894687916)
               BPARA = 1.0129458
               CPARA = -0.604615023/(1+9.5779826*exp(-1*2.2080809*ld))
               cd = APARA*(BPARA**mach)*(mach**(CPARA))
            elif mach >= 1.05 and mach <= 1.2:
               MPARA = 1/(-0.156531249-0.35165656*ld)
               BPARA = 1/(0.10668068+0.2160142549*ld)
               cd = MPARA*mach+BPARA
            elif mach <1.05:
               cd = -0.075*ld + 0.275

        if (self.Nosecone == 'parabolic'): # assume that both nosecone shapes give same shape
            if mach >= 1.2: # All values taken from fits to HyperCFD
               APARA = 0.2433566382*exp(ld**-7.1807129)
               BPARA = 1.009709
               CPARA = -0.567521484056/(1+5.59560038938568*exp(-1*2.23635526782648*ld))
               cd = APARA*(BPARA**mach)*(mach**(CPARA))
            elif mach >= 1.05 and mach <= 1.2:
               MPARA = 1/(-0.1595385088-0.41826608398336*ld)
               BPARA = 1/(0.123489761128+0.244747303231711*ld)
               cd = MPARA*mach+BPARA
            elif mach <1.05:
               cd = -0.025*ld + 0.125


        if (self.Nosecone == 'sears-haack'): # assume that both nosecone shapes give same shape
            if mach >= 1.2: # All values taken from fits to HyperCFD
               APARA = 0.243884345*exp(ld**-0.80690309)
               BPARA = 1.0047095
               CPARA = -0.60330669/(1+14.6196741884*exp(-1*3.27801239521*ld))
               cd = APARA*(BPARA**mach)*(mach**(CPARA))
            elif mach >= 1.05 and mach <= 1.2:
               MPARA = 1/(-0.111417758-0.436291862*ld)
               BPARA = 1/(0.090907066+0.26210278132*ld)
               cd = MPARA*mach+BPARA
            elif mach <1.05:
               cd = -0.05*ld + 0.25

		
        cd = cd +cdfins
        #print 'cd,cdfins',cd,cdfins,'    mach,h',mach,h/1000.    
        return cd
        

    def to_radians(self,degree):
        return degree * pi/180
        
if __name__ == "__main__":
    print "the simulation object"
    print "using simple text interface, minimum energy trajectory"
    print ""
    sim = Simulation(None) #this simulation object has no parent
    sim.numstages = int(raw_input("Number of stages: "))
    sim.numfins = int(raw_input("Number of fins: "))
    for i in range(1,sim.numstages+1):
        sim.fuelmass.append(float(raw_input("Fuel mass: ")))
        drymass = (float(raw_input("Dry mass: ")))
        sim.m0.append(drymass + sim.fuelmass[i])
        sim.fuelfraction.append(sim.fuelmass[i]/sim.m0[i])
        sim.Isp0.append(float(raw_input("Isp: ")))
        sim.thrust0.append(float(raw_input("Thrust (kg f): "))*9.81)
        sim.dMdt.append(float(sim.thrust0[i]/(sim.Isp0[i]*9.81)))
        sim.burntime.append(float(raw_input("Burntime (sec): ")))
        sim.payload = float(raw_input("Payload (kg): "))
        sim.missilediam = float(raw_input("Missile Diameter (m): "))
        sim.nozzlearea = float(raw_input("Nozzle Area (m^2): "))
        sim.LdivD = float(raw_input("Nosecone: Length/Diam (dimensionless): ")) 
        sim.rvdiam = float(raw_input("Re-entry Diameter (m): "))
        sim.est_range = float(raw_input("Est range (km): "))*1000
        sim.burnout_angle = float(raw_input("Burnout Angle (deg): "))*1


        # ADDED MAY 16
        sim.FinHeight = float(raw_input("Fin Height [m]: "))*1
        sim.SweepAngle = float(raw_input("Sweep Angle [deg]: "))*1
        sim.FinTipChordLength = float(raw_input("Fin Tip Chord Length [m]: "))*1
        sim.FinRootChord = float(raw_input("Fin Root Chord [m]: "))*1
        sim.FinThickness = float(raw_input("Fin Thickness [m]: "))*1
        #

    print '\n'
    sim.trajectory = "Minimum Energy"
    results = sim.integrate(sim.trajectory)
    print '\n'
    
    path = 'data.txt'
    outfile = open(path,'w')
    for i in range(1,sim.numstages+1):
                outfile.write("STAGE %i Parameters:\n" % i)
                outfile.write("Fuel mass (kg): " + str(sim.fuelmass[i]) + '\n')
                outfile.write("Dry mass (kg): " + str(sim.m0[i] - sim.fuelmass[i]) + '\n')
                outfile.write("Fuel fract: " + str(sim.fuelfraction[i]) + '\n')
                outfile.write("Isp @ SL: " + str(sim.Isp0[i]) + '\n')
                outfile.write("Burn time (sec): " + str(sim.burntime[i]) + '\n')
                outfile.write("Thrust (N): " + str(sim.thrust0[i]) + '\n')
                outfile.write("dM/dt: " + str(sim.dMdt[i]) + '\n')
                
    outfile.write("\nTIME,HEIGHT,VELOCITY,MASS,THRUST,DRAG,CD,GAMMA,RANGE\n")
    
    flat = zip(results['Time'],
                results['Height'],
                results['Velocity'],
                results['Mass'],
                results['Thrust'],
                results['Drag'],
                results['Cd'],
                results['Gamma'],
                results['Range'])
    for i in range(1,len(flat)):
        for n in range(0,len(flat[i])):
            outfile.write('%.3f' % flat[i][n])
            outfile.write(',')
        outfile.write('\n')
    print "Data written to '%s'" % path
    self.outfile.close()

else:
    import wx
    #using gui
