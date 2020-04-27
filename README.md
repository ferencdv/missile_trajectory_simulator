# missile_trajectory_simulator
This code is directly taken from Josh Levinger's Missile Trajectory simulator code. 
Modifications are by me only and will be to improve its utility for our applications. 
All credit should go to Josh Levinger of GlobalSecurity.org and Union of Concerned Scientists 
nuclear expert Dr. Dave Wright on which this code is based. You can find Josh's webpage and 
code at: https://github.com/jlev/ballistic-missile-range. Also his Readme page contains 
all relevant information as to the licensing of this code.

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

April 27 2020: The presets file has also been changed to add the Iranian Qasam missile with and without coasting. Coasting is included as an empty stage with zero thrust. Specifically if ISP=1, and Mass Fuel = 0.01 kg, then in the burn time is t seconds then the thrust = (MF/t) = 0.01/42 = 0.000238 kgf. Similarly if t = 200 s then thrust = MF/t = 0.01/200 = 0.00005 kgf.

April 27 2020: Note that the uptodate files are: gui-11.py, presets.gui11.py, simss.py
