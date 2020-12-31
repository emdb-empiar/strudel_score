# strudel_score
A tool for visualising the STRUDEL map-model validation results. Runs as a ChimeraX plugin.


Installation:

1.    Install ChimeraX
2.    Unpack strudel_score.tar
3.    From command line go to strudel_score folder
4.    On mac and linux type: 

     make install 
     
     
   (If it gives an error please check the path to Chimera in Makefile, line 17)
   
   On windows: 
  
     make_win.bat app-install

Usage:

1.    Open ChimeraX
2.    Open Strudel Score: Tools/General/Strudel Score
3.    Set the motifs library path. 
      
      Download and unpack:
      
      https://drive.google.com/drive/folders/1G0-XQOxCUHkcBCA2CQPSfkPLiPnEA4fw?usp=sharing
      
      In the Strudel Score window: File/Set Motif Library
	  and select the unpacked folder with strudel_libs
	  
4. Download validation results examples:

    https://drive.google.com/drive/folders/1upXLVZf7kout0mcgt5iqabonvRM9FWih?usp=sharing
   
   In the Strudel Score window: File/Open Project and select an unpacked validation results folder
	  