# MISSILE TRAJECTORY SIMULATOR
This code is directly taken from Josh Levinger's Missile Trajectory simulator code. 
Modifications are by me only and will be to improve its utility for our applications. 
All credit should go to Josh Levinger of GlobalSecurity.org and Union of Concerned Scientists 
nuclear expert Dr. Dave Wright on which this code is based. You can find Josh's webpage and 
code at: https://github.com/jlev/ballistic-missile-range. Also his Readme page contains 
all relevant information as to the licensing of this code.

## Inclusion of Drag
May 14

Note that this is a fit of total drag which includes the body drag plus the fricton drag 
and base drag. It does not include the fin drag. Note that this drag also does not include 
the drag due to various length bodies because it does not seem to be a large effect. 
The HyperCFD code is developed by John Cipolla and can be obtained from his website at:
http://www.aerorocket.com/hyperx.html. I have used his results essentially as a look up
table.
    
### Fitting Procedure

The drag coefficient CD functions for various L/D (length of nosecone/ max diameter of 
nosecone) and for various mach speeds were produced using HyperCFD code. Note that this 
is a fit of total drag which includes the body drag plus the fricton drag and base drag. 
It does not include the fin drag. Note that this drag also does not include 
the drag due to various length bodies because it does not seem to be a large effect. 
These curves can be divided into three segments from M=0 to M=1.05 [1], from M=1.05 
to M=1.2 [2] and from M=1.2 and above. In each of these sections the curves vary as a 
function of M and of L/D ratio. Segment [1] is the easiest and is interpolated linearly
from L/D=1 to L/D=3. The next segment [2] is a very steep fall off to M=1.2 and 
is also intepreted linearly from M=1.05 to M=1.2 (essentially a linear interpolation
beteen these two points. 

![HyperCFD CD as a function of Mach](https://github.com/ferencdv/missile_trajectory_simulator/blob/master/Capture.JPG)

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

---------------------------------------------------------------------------------------------

April 23 2020: Please use new files. gui-7-new.py which also requires the simss.py files.
April 24 2020: Still have error in valueerror popup window in gui-11.py (this is the latest file and it depends on simss.py)
           dlg = wx.MessageDialog(self,"Please make sure all fields are filled in.","Entry error",wx.OK | wx.ICON_INFORMATION)
I can prevent it from giving this eror by changing:
            sim.numstages = int(self.StageChoiceBox.GetSelection()+1) to sim.numstages = int(self.StageChoiceBox.GetSelection()+0)
 But the problem is that this does not properly populate the results in the 'results panel' the next panel. Also, the results are not correct, so I think it does not incluse a stage in the calculations.
 
April 27 2020: There was a major bug in the original program in that when a new stage would start to burn the mass was set back to the original mass of the missile. This was fixed in simss.py with the following code:

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
                
 This is added before the END OF BIG LOOP. This code redefines the mass at each stage to account for a loss of a stage. Note that this this code happens at the end of the stage integration.

April 27 2020: The presets file has also been changed to add the --- missile with and without coasting. Coasting is included as an empty stage with zero thrust. Specifically if ISP=1, and Mass Fuel = 0.01 kg, then in the burn time is t seconds then the thrust = (MF/t) = 0.01/42 = 0.000238 kgf. Similarly if t = 200 s then thrust = MF/t = 0.01/200 = 0.00005 kgf.

April 27 2020: Note that the uptodate files to use/download are: gui-11.py, presets.gui11.py, simss.py
