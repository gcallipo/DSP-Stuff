# DSP stuff
This project define a space on which to share components, algorithms and exercises concerning 
the world of digital signal analysis, this includes SDR techniques, Audio and Video manipulation, 
controls etc.

#### Included Projects:
	
	Audio Filter Convolutional F32: 
	
	This is a specific audio signal filtering technique implemented 
	in the frequency domain.
	
	This software is based on the  AudioFilterConvolution routine 
	Written by Brian Millier on Mar 2017, more infos here:
	https://circuitcellar.com/research-design-hub/fancy-filtering-with-the-teensy-3-6/
	
	and modified by Giuseppe Callipo-ik8yfw (myself)
	
	Modification:
	
		* Class refactoring, change some methods visibility;
		* Filter coefficients calculation included into class;
		* Change the class for running in both with F32  OpenAudio_ArduinoLibrary for Teensy;	
		* Added initFilter method for single anf fast initialization and on the fly reinititializzation; 
		* Optimize it to use as output audio filter on SDR receiver.
		
		
#### NOTE:

All software included in the projects are intended as opensources.
(Updated on 02.05.2021)
