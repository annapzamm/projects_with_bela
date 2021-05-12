/*
 ____  _____ _        _    
| __ )| ____| |      / \   
|  _ \|  _| | |     / _ \  
| |_) | |___| |___ / ___ \ 
|____/|_____|_____/_/   \_\
The platform for ultra-low latency audio and sensor processing
http://bela.io
A project of the Augmented Instruments Laboratory within the
Centre for Digital Music at Queen Mary University of London.
http://www.eecs.qmul.ac.uk/~andrewm
(c) 2016 Augmented Instruments Laboratory: Andrew McPherson,
	Astrid Bin, Liam Donovan, Christian Heinrichs, Robert Jack,
	Giulio Moro, Laurel Pardue, Victor Zappi. All rights reserved.
The Bela software is distributed under the GNU Lesser General Public License
(LGPL 3.0), available here: https://www.gnu.org/licenses/lgpl-3.0.txt
*/

#include <Bela.h>
#include <cmath>
//#include <libraries/Scope/Scope.h>
#include <libraries/WriteFile/WriteFile.h> 
#include <libraries/Trill/Trill.h>
#include <vector>
#include <iostream>
#include <string>
#include <string.h>



std::vector<Trill*> touchSensor; // Trill object declaration
#define NUM_TOUCH 5 // Number of touches on Trill sensor
#define NUM_SENSORS 2 // Number of sensors. 

// Location of touches on Trill Ring
double gTouchLocation[NUM_SENSORS][NUM_TOUCH] = {{ 0.0, 0.0, 0.0, 0.0, 0.0 }, 
	{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
	
// Size of touches on Trill Ring
float gTouchSize[NUM_SENSORS][NUM_TOUCH] = {{ 0.0, 0.0, 0.0, 0.0, 0.0 }, 
	{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
	
// Number of active touches
unsigned int gNumActiveTouches[NUM_SENSORS] = {0, 0};
double read_out[NUM_SENSORS] = {0.0, 0.0};
double past_read[NUM_SENSORS] = {0.0, 0.0};
double past_read_out[NUM_SENSORS] = {0.0, 0.0}; 
float active_touch[NUM_SENSORS] = {0.0, 0.0};
float num_unique_trill_locs = 3584.0; 

// tone characteristics
float gPhase[NUM_SENSORS] = {0.0, 0.0}; // carrier phase
float gInverseSampleRate; // sampling period. 
float gAmplitude= 0.01; // amplitude.
float centerFrequency[NUM_SENSORS] = {0.0, 0.0}; // peak frequency (Hz) at 0/360 
float freqRange_aroundCenter = 500.0;//100.0; // range around center frequency (180 degrees = centerFrequency + freqRange_aroundCenter, Hz)
double frequency[NUM_SENSORS] = {0.0, 0.0};
double read_out_degrees[NUM_SENSORS] = {0.0, 0.0}; 
double y_coordinate_now[NUM_SENSORS] = {0.0, 0.0}; 
double pastRead[NUM_SENSORS] = {0.0, 0.0};
double newRead[NUM_SENSORS] = {0.0, 0.0}; 	


// sleep time, etc.
unsigned int sample_num = 0;
unsigned int gTaskSleepTime  = 2500; 

// Browser-based oscilloscope
//Scope gScope;

//output file.
WriteFile freq_mod_tone_bela; 
WriteFile tone_vals;
//ofstream outdata; // outdata is like cin

// timstamp
float timestamp_ms = 0.0;


//trial duration. 
#define TRIALDUR_MS 60000 

// input arguments (passed from bash)
std::string name; // output data filename
std::string exp_version; // exp version (discrete / continuous)
std::string aud_feedback_cond; //auditory feedback condition. 
std::string octave_order; //octave order. 

// metronome.
#include <math.h>

// Oscillator variables
float gPhaseMet = 0;			// Current phase
float gFrequencyMet = 1000;	// Frequency in Hz

// metronome variables.
float gAmplitudeMet = 0;   
float gAmplitudeMet_peak = 0.01;
float gEnvelopeScaler = 0.997;
float gMetronomeCounter = 0;
float gMetronomeInterval = 0;
float gMetClickNum = 0; 
float gMetTrillLocCounter = 0; 
float touch_interval_met = 0.0; 
double met_trill_pos_orig = 0.0; 
double met_freq_trill[NUM_SENSORS] = {0.0, 0.0}; 
float gPhaseMetContinuous[2] = {0.0, 0.0}; 
float ramp_increment = 0.001;
float ramp_counter = 0.0;
float ramp_duration_samples = 441; // at 44.1 kHz sampling rate.


void loop(void*)
{
	
	std::vector<float> logs(500000);
	unsigned int writeFileSize = logs.size() + 5000;
	freq_mod_tone_bela.setBufferSize(writeFileSize);
	unsigned int logIdx = 0;

	
	while(!Bela_stopRequested())
	{

		// Read locations from Trill sensor
		unsigned int q = touchSensor.size() < NUM_TOUCH ? touchSensor.size() : NUM_TOUCH;
		float timestamp_now = timestamp_ms; // same timestamp for both sensors.

		for(unsigned int n = 0; n < q; ++n){
			
			float data_out_now[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			data_out_now[0] = timestamp_now;
			data_out_now[1] = float(n);
			

			// read current sensor & number of active touches.
			touchSensor[n]->readI2C();
			gNumActiveTouches[n] = touchSensor[n]->getNumTouches();

			// loop thru active touches, get size & location.
			for(unsigned int i = 0; i < gNumActiveTouches[n]; i++) {
				gTouchLocation[n][i] = touchSensor[n]->touchLocation(i);
				gTouchSize[n][i] = touchSensor[n]->touchSize(i);
			}
			// For all inactive touches, set location and size to 0
			for(unsigned int i = gNumActiveTouches[n]; i < NUM_TOUCH; i++) {
				gTouchLocation[n][i] = 0.0;
				gTouchSize[n][i] = 0.0;
			}
			// if there are active touches.
			if(touchSensor[n]->getNumTouches())
			{
				past_read[n] = read_out[n]; // store past read.
				read_out[n] = gTouchLocation[n][0]; // get touch location for first touch (of 5).
				active_touch[n]=1.0; // variable to tell audio loop whether there is active touch.

				// Subtract .25 from touch location & wrap between 0 & 1.
				// (this places 0 degrees on x-axis, i.e. at original 270 degrees in Trill schematic)
				// change made on Feb 26 2021
				read_out[n] = read_out[n] + 0.25; 
				if (read_out[n] > 1.0 ){
					read_out[n] = read_out[n] - 1.0; 
				}

				// now convert position to value between 0 & 360.
				read_out_degrees[n] = read_out[n] * 360.0;

				// now compute sine of that value to get y-coordinate.
				// if using bar sensor, don't use y-value.
			    y_coordinate_now[n] = sin(read_out_degrees[n]*(M_PI/180.0));// if bar sensor, motion is along x-axis.

				// determine frequency.
				frequency[n] = map(y_coordinate_now[n], -1.0, 1.0, centerFrequency[n]-freqRange_aroundCenter, centerFrequency[n]+freqRange_aroundCenter);
				
				//rt_printf("frequency range for sensor %d is %f-%f\n", n, centerFrequency[n]-freqRange_aroundCenter, centerFrequency[n]+freqRange_aroundCenter);
				rt_printf("read out sensor %d is %f\n", n, frequency[n]);
			
				// log sensor data IF new touch location.
				if(read_out[n]!=past_read[n]){
					data_out_now[0] =	timestamp_ms;
					data_out_now[1] = float(n);
					data_out_now[2] = read_out[n];
					data_out_now[3] = y_coordinate_now[n];
					data_out_now[4] = frequency[n];
					data_out_now[5] = 1; 
				}
			} else {
				active_touch[n]=0.0;
			}
			
			 
		
		for(unsigned int n1 = 0; n1 < 10; ++n1)
				{
					logs[logIdx++] = data_out_now[n1]; // fill an array in memory
					
					if(logs.size() >= writeFileSize){
						rt_printf("Data has exceeded output file memory!!!!!!!!!!!!");
					}

						//if(logs.size() == logIdx)
						if(timestamp_ms >= TRIALDUR_MS )
						{
							// when the array is full, dump it to disk
							// given how we have setBufferSize() above, WriteFile's buffer
							// should have enough space for all of logs
							printf("Dumping\n");
							freq_mod_tone_bela.log(logs.data(), logIdx);
							logIdx = 0;
							//for(unsigned int l1 = 0; l1 < logIdx; ++l1){
							//outdata << logs[l1] << endl;
							//}
							// now we need to wait until the WriteFile thread has written everything to disk
							// unfortunately it will always keep up to 4096 samples inside its buffer
							while(freq_mod_tone_bela.getBufferStatus() < (writeFileSize - logs.size()) / writeFileSize)
							{
								printf("waiting %f\n", freq_mod_tone_bela.getBufferStatus());
								usleep(10000);
							}
							printf("stopping\n");
							Bela_requestStop();
						}
					}
		}

	usleep(gTaskSleepTime);
	}
	// send WriteFile any leftover stuff
	freq_mod_tone_bela.log(logs.data(), logIdx);
	//outdata.close();
}


bool setup(BelaContext *context, void *userData)
{
	
	// get input arguments from bash.
	std::cin >> name >> exp_version >> aud_feedback_cond >> octave_order;
	
	//now summarize in terminal. 
    std::cout << "Filename is: "; //filename.
    std::cout << name << std::endl;
    
    std::cout << "Exp version is: "; //exp version.
    std::cout << exp_version << std::endl;
    
    std::cout << "Auditory feedback condition is: ";
    std::cout << aud_feedback_cond << std::endl; 
    
    
    std::cout << "Octave order is : ";
    std::cout << octave_order << std::endl; 
    
    // if octave order = 0, then subject with lower number is in lower octave. 
    if(std::string(octave_order) == "0"){ 
    	centerFrequency[0] = 1000.0;
    	centerFrequency[1] = 1500.0;
    } else {
    	rt_printf("NOW");
    	centerFrequency[0] = 1500.0;
    	centerFrequency[1] = 1000.0;
    }
    
    char const *name_out = name.data(); 
    //char const *exp_version_out = exp_version.data(); 
    char filename_cat_out[strlen("/root/Bela/projects/anzjal_exp1_trill_continuous/") + strlen(name_out) + 1];
    strcpy(filename_cat_out, "/root/Bela/projects/anzjal_exp1_trill_continuous/");
    strcat(filename_cat_out, name_out);

	// output file.
	freq_mod_tone_bela.setup(filename_cat_out);
	freq_mod_tone_bela.setFormat("%f %f %f %f %f %f %f %f %f %f\n");
	freq_mod_tone_bela.setFileType(kBinary); 
		
	// output file - stimulus frequency mapping.
	//char filename_cat_out_stim[strlen("/root/Bela/projects/anzjal_exp1_trill_continuous/") + strlen("stimulus_mapping.m") + 1];
	//strcpy(filename_cat_out_stim, "/root/Bela/projects/anzjal_exp1_trill_continuous/");
    //tone_vals.setup(filename_cat_out_stim);
    //tone_vals.setFormat("%f %f %f %f %f %f\n");
    //tone_vals.setFileType(kText);
    
	gInverseSampleRate = 1.0 / context->audioSampleRate;
	
	// set up Trill sensors. 
	unsigned int num_sensors=0;
	unsigned int i2cBus = 1;
	for(uint8_t addr = 0x20; addr <= 0x50; ++addr)
	{
		Trill::Device device = Trill::probe(i2cBus, addr);
		if(Trill::NONE != device && Trill::CRAFT != device)
		{
			num_sensors++;
			touchSensor.push_back(new Trill(i2cBus, device, addr));
			printf("%#4x (%3d) | %s\n", addr, addr, Trill::getNameFromDevice(device).c_str());
			touchSensor.back()->setScanSettings(0,10); // check this. 
			touchSensor.back()->setAutoScanInterval(1);
			touchSensor.back()->printDetails();
		}
	}

	//touchSensor.printDetails();
	
	// Set and schedule auxiliary task for reading sensor data from the I2C bus
	Bela_runAuxiliaryTask(loop, 90); //high priority
	
	// Set up the oscilloscope
	//gScope.setup(NUM_SENSORS, context->audioSampleRate);
	
	// metronome.
	// Calculate the metronome interval in samples based on 75 bpm
	float bpm = 46.14258;
	gMetronomeInterval = 60.0 * context->audioSampleRate / bpm;
	
	// number of samples between updating touch locations for 1 cycle around Trill at metronome rate. 
	touch_interval_met = (float)gMetronomeInterval / (float)num_unique_trill_locs; 


	
	
	return true;
}

void render(BelaContext *context, void *userData)
{
	// define output amplitude (default zero unless active touch).
	float out_amp[2] = {0.0, 0.0};

	// loop thru samples in audio block. 
	for(unsigned int n = 0; n < context->audioFrames; n++) {
    
		    float outMet[2] = {0.0, 0.0}; 
		    float met_trill_pos = 0.0; 
		    
		   // now on each sample since turn sound on, increment value.
			 if(gMetClickNum  < 7){
    			if(gAmplitudeMet < gAmplitudeMet_peak){
        			gAmplitudeMet = gAmplitudeMet + ramp_increment;
    			}
    		} else if (gMetClickNum == 7 && ((gMetronomeInterval-gMetronomeCounter)==ramp_duration_samples)){
	    		if(gAmplitudeMet > 0){
	        		gAmplitudeMet = gAmplitudeMet - ramp_increment;
	    		}
    		}
    		
	        // if time to update metronome finger position & corresponding frequencies, then do it. 
	        if(gMetTrillLocCounter==0){
	        	
				  //update trill position in metronome cycle.
			      if(gMetronomeCounter > 0){
			        met_trill_pos_orig = met_trill_pos_orig + (double)(1.0/num_unique_trill_locs);
			      } else {
			        met_trill_pos_orig = 0;
			      }
			      
			      
			      //Now store original value.
			      met_trill_pos = met_trill_pos_orig;
			      
			      // Now convert position to value between 0 & 360.
			      double met_trill_pos_degrees = met_trill_pos * 360.0;
			      
			      // Now compute y-coordinate of position.
			      double y_coordinate_met = sinf(met_trill_pos_degrees * (M_PI/180.0));
			      
			      // determine corresponding frequency.
			      // map y-coordinates -.5:.5 linearly onto min frequency:max frequency in hz.
			      met_freq_trill[0] = map(y_coordinate_met, -1.0, 1.0, centerFrequency[0]-freqRange_aroundCenter, centerFrequency[0]+freqRange_aroundCenter);
			      met_freq_trill[1] = map(y_coordinate_met, -1.0, 1.0, centerFrequency[1]-freqRange_aroundCenter, centerFrequency[1]+freqRange_aroundCenter);
			      
			      
			      //float output_mapping[6] = {(float)timestamp_ms, (float)met_trill_pos_orig, (float)met_trill_pos_degrees, (float)y_coordinate_met, (float)met_freq_trill[0], (float)met_freq_trill[1]};
			      //if(gAmplitudeMet > 0){
			      // tone_vals.log(output_mapping, 6);
			      //}
			   
	        }
	        
	        
				
		
	    	// Calculate a sample of the sine wave, and scale by the envelope
	    	if(std::string(exp_version)=="discrete"){
	    		rt_printf("true");
				gPhaseMet += 2.0 * M_PI * gFrequencyMet / context->audioSampleRate;
				if(gPhaseMet >= 2.0 * M_PI)
					gPhaseMet -= 2.0 * M_PI;
				outMet[0] = gAmplitudeMet * sin(gPhaseMet);
				outMet[1] = gAmplitudeMet * sin(gPhaseMet);
				
		    	} else {
		    		
		    		for(unsigned int t1 = 0; t1 < NUM_SENSORS; ++t1){
		    			gPhaseMetContinuous[t1] += 2.0 * M_PI * met_freq_trill[t1] / context->audioSampleRate;
		    			if(gPhaseMetContinuous[t1] >= 2.0 * M_PI)
							gPhaseMetContinuous[t1] -= 2.0 * M_PI;
						outMet[t1] = gAmplitudeMet * sin(gPhaseMetContinuous[t1]);
						
						
		    		}
		    	}
	    	
		
		    // sample count.
			sample_num++;
		
			// define array for sensor audio output.
			float out_freqs[NUM_SENSORS] = {0.0, 0.0};
			
			// only need to record finger location / frequency at first sample in audio block, as does not change during block. 
			timestamp_ms = float(sample_num) * (1000/(context->audioSampleRate)); 
		
			// loop thru Trill sensors & compute phase. 
			for(unsigned int t = 0; t < NUM_SENSORS; ++t) {
						
					gPhase[t] += 2.0f * (float)M_PI * frequency[t] * gInverseSampleRate;
					if( gPhase[t] > M_PI ){
						gPhase[t] -= 2.0f * (float)M_PI;
					}
					
					// set output amplitude ON if there is active touch AND if trial duration is not exceeded. 
					if(timestamp_ms < TRIALDUR_MS){ 
						if(active_touch[t] > 0){
							//if (gMetClickNum >= 8){
							 out_amp[t] = 0.01; 
							//}
						}
					}
					
					// now store output in array. 
					out_freqs[t] =  out_amp[t] * sinf(gPhase[t]);
				
					
				
			} // end looping through sensors. 
			
			
			unsigned int sensor_num_audio_self;
			unsigned int sensor_num_audio_partner;
			
			for(unsigned int channel = 0; channel < context->audioOutChannels; channel++) {
				
				float out = 0.0;
				
				// send out sensor 1 data on channel 1, sensor 2 data on channel 2.
				if(channel==0){
					sensor_num_audio_self=1;
					sensor_num_audio_partner=0;
				} else {
					sensor_num_audio_self=0;
					sensor_num_audio_partner=1;
				}
				
				// if metronome signal is over.
				if (gMetClickNum >= 8){ 
				
					// define audio output by condition. 
					// uncoupled = 0 > 0; 1 > 1
					// unidirec1 = 0 > 1; 1 > 1 
					// unidirec2 = 0 > 0; 1 > 0
					// bidirect = 0 > 1; 1 > 0
					// metronome pitch always corresponds to pitch of whatever channel subject is hearing. 
					
					
					if(std::string(aud_feedback_cond) == "uncoupled" | std::string(aud_feedback_cond) == "exploration" | std::string(aud_feedback_cond) == "practice"){
						
						out+=out_freqs[sensor_num_audio_self];

					} else if (std::string(aud_feedback_cond) == "uni1"){
						if(channel==0){
							
							out+=out_freqs[sensor_num_audio_partner];

						} else {
							
							out+=out_freqs[sensor_num_audio_self];

						}
						
					} else if (std::string(aud_feedback_cond) == "uni2"){
						if(channel==1){
							
							out+=out_freqs[sensor_num_audio_partner];

						} else {
							
							out+=out_freqs[sensor_num_audio_self];

						}
					} else if (std::string(aud_feedback_cond) == "bidirectional"){
						
						out+=out_freqs[sensor_num_audio_partner];

					}
			    // else if metronome is not over, then only play metronome UNLESS exploration phase (in which case just allow subjects to move around sensor).
				} else {
					
					if(std::string(aud_feedback_cond) == "uncoupled" | std::string(aud_feedback_cond) == "practice"){
						
						out+=outMet[sensor_num_audio_self];
						
					} else if (std::string(aud_feedback_cond) == "uni1"){
						if(channel==0){
							
							out+=outMet[sensor_num_audio_partner];
							
						} else {
							
							out+=outMet[sensor_num_audio_self];
							
						}
						
					} else if (std::string(aud_feedback_cond) == "uni2"){
						if(channel==1){
							
							out+=outMet[sensor_num_audio_partner];
							
						} else {
							
							out+=outMet[sensor_num_audio_self]; 
							
						}
					} else if (std::string(aud_feedback_cond) == "bidirectional"){
						
						out+=outMet[sensor_num_audio_partner];
						
					} else if (std::string(aud_feedback_cond) == "exploration"){
						out+=out_freqs[sensor_num_audio_self];

					}
					
				}
				
			   
			   // now write out audio. 
			   audioWrite(context, n, channel, out);
			    	
			}
			
			
		//gScope.log(out_freqs[0], out_freqs[1]);
		 //////////metronome.
		    // Increment counter and check if it has reached the interval
			if(++gMetronomeCounter >= gMetronomeInterval) {
				// Metronome tick elapsed; reset counter and envelope
				gMetronomeCounter = 0;
				met_trill_pos_orig = 0.0; 
				gMetTrillLocCounter = 0; 
				//gAmplitudeMet =  0.8;
				++gMetClickNum;
			}
		
		    //if(gMetClickNum >= 8){ 
		    //	gAmplitudeMet = 0; 
		    //} 
	
	        // reset every 6 cycles. 
	        if(++gMetTrillLocCounter >= touch_interval_met){
	        	gMetTrillLocCounter = 0;
	        }
	        
		
	}

}

void cleanup(BelaContext *context, void *userData)
{

}