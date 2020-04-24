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
