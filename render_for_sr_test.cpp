
#include <Bela.h>
#include <cmath>
#include <libraries/Scope/Scope.h>
#include <libraries/WriteFile/WriteFile.h> 
#include <libraries/Trill/Trill.h>
#include <vector>

std::vector<Trill*> touchSensor; // Trill object declaration
#define NUM_TOUCH 5 // Number of touches on Trill sensor
#define NUM_SENSORS 2 // Number of sensors. 
// Location of touch on Trill Ring
float gTouchLocationCycle[NUM_SENSORS] = {0.0, 0.0};
// Size of touch on Trill Ring
float gTouchSizeCycle[NUM_SENSORS] = {0.0, 0.0};

// Location of touches on Trill Ring
float gTouchLocation[NUM_SENSORS][NUM_TOUCH] = {{ 0.0, 0.0, 0.0, 0.0, 0.0 }, 
	{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
// Size of touches on Trill Ring
float gTouchSize[NUM_SENSORS][NUM_TOUCH] = {{ 0.0, 0.0, 0.0, 0.0, 0.0 }, 
	{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
// Number of active touches
unsigned int gNumActiveTouches[NUM_SENSORS] = {0, 0};
//unsigned int gTaskSleepTime = 12000; // microseconds
float read_out[NUM_SENSORS] = {0.0, 0.0};
float past_read[NUM_SENSORS] = {0.0, 0.0};


float gPhase[NUM_SENSORS] = {0.0, 0.0}; // carrier phase
float gInverseSampleRate; // sampling period. 
float gAmplitude=.6; // amplitude.
float peakFrequency[NUM_SENSORS] = {1000.0, 1000.0}; // peak frequency (Hz) at 0/360 
float freqRange_aroundPeak = 250;//100.0; // range around center frequency (180 degrees = peakFrequency + freqRange_aroundPeak, Hz)
float past_read_out[NUM_SENSORS] = {0.0, 0.0}; 
float active_touch[NUM_SENSORS] = {0.0, 0.0};
unsigned int sample_num = 0;
unsigned int gTaskSleepTime  = 2500; 

// Browser-based oscilloscope
Scope gScope;

//output file.
WriteFile freq_mod_tone_bela; 

float frequency[NUM_SENSORS] = {0.0, 0.0};
float read_out_degrees[NUM_SENSORS] = {0.0, 0.0}; 
float y_coordinate_now[NUM_SENSORS] = {0.0, 0.0}; 
float timestamp_ms = 0.0;
float pastRead[NUM_SENSORS] = {0.0, 0.0};
float newRead[NUM_SENSORS] = {0.0, 0.0}; 	



void loop(void*)
{

	while(!Bela_stopRequested())
	{
		
		  float data_out_now[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		  data_out_now[0] = { timestamp_ms };
		  freq_mod_tone_bela.log(data_out_now,10);
		     
			// Read locations from Trill sensor
			for(unsigned int n = 0; n < NUM_SENSORS; ++n){
				
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
					
					// now convert position to value between 0 & 360. 
					read_out_degrees[n] = read_out[n] * 360.0; 
					
					// now compute sine of that value to get y-coordinate.
					// if using bar sensor, don't use y-value. 
					if(n==1){
						y_coordinate_now[n] = sinf(read_out_degrees[n]*M_PI/180.0) * 0.5; //NECESSARY WITH RING. 
					} else {
						y_coordinate_now[n] = cosf(read_out_degrees[n]*M_PI/180.0) * 0.5;// if bar sensor, motion is along x-axis. 
					}
	
					// determine frequency. 
					frequency[n] = map(y_coordinate_now[n], -0.5, 0.5, peakFrequency[n]-freqRange_aroundPeak, peakFrequency[n]+freqRange_aroundPeak);
					
					// log sensor data IF new touch location.
					//if(read_out[n]!=past_read[n]){
						//float data_out_now[5] =	{ timestamp_ms, float(n), read_out[n], y_coordinate_now[n], frequency[n] };
						//freq_mod_tone_bela.log(data_out_now,5);
						//float tmp_az[10] = { timestamp_ms, float(n), read_out[n], y_coordinate_now[n], frequency[n], gTouchLocation[n][0],  gTouchLocation[n][1] ,  gTouchLocation[n][2] ,  gTouchLocation[n][3] ,  gTouchLocation[n][4]  };
						//for(unsigned int tmp_i = 0; tmp_i < 10; tmp_i++){
							//data_out_now[tmp_i] = tmp_az[tmp_i];
						//}
						//freq_mod_tone_bela.log(data_out_now,10);
						//freq_mod_tone_bela.log(data_out_now,10);
					//}
				} else {
					active_touch[n]=0.0;
					//if(n==1){
				     // data_out_now[0] = { timestamp_ms };
					//}
				}
				
		}
	
		usleep(gTaskSleepTime);
	}
}


bool setup(BelaContext *context, void *userData)
{
	
	
	// output file.
	freq_mod_tone_bela.setup("freq_mod_tone_bela.bin");
	//freq_mod_tone_bela.setFormat("%f %f %f %f %f\n");
	freq_mod_tone_bela.setFormat("%f %f %f %f %f %f %f %f %f %f\n");
	freq_mod_tone_bela.setFileType(kBinary); 
	
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
	Bela_runAuxiliaryTask(loop);
	
	// Set up the oscilloscope
	gScope.setup(NUM_SENSORS, context->audioSampleRate);
	
	// generate sine tone representing modulation frequency range. 
	
	

	
	return true;
}

void render(BelaContext *context, void *userData)
{
	// define output amplitude (default zero unless active touch).
	float out_amp = 0.0;

	// loop thru samples in audio block. 
	for(unsigned int n = 0; n < context->audioFrames; n++) {
		
		    // sample count.
			sample_num++;

			// define array for audio output.
			float out_freqs[NUM_SENSORS] = {0.0, 0.0};
			
			// only need to record finger location / frequency at first sample in audio block, as does not change during block. 
			timestamp_ms = float(sample_num) * (1000/(context->audioSampleRate)); 
		
			// loop thru Trill sensors & compute phase. 
			for(unsigned int t = 0; t < NUM_SENSORS; ++t) {
						
					gPhase[t] += 2.0f * (float)M_PI * frequency[t] * gInverseSampleRate;
					if( gPhase[t] > M_PI ){
						gPhase[t] -= 2.0f * (float)M_PI;
					}
					
					// set output amplitude ON if there is active touch. 
					if(active_touch[t] > 0){
						 out_amp = 0.8; 
					}
					
					// now store output in array. 
					out_freqs[t] =  out_amp * sinf(gPhase[t]);
				
			} // end looping through sensors. 

			// loop thru audio channels. 
			for(unsigned int channel = 0; channel < context->audioOutChannels; channel++) {
				
				float out = 0;
				float stream_now=0;
				
				for(int i=0;i<NUM_SENSORS;i++){
					if(i==0){
						stream_now=out_freqs[i];
					} else if (i==1){
						stream_now=out_freqs[i];
					}
				    out += stream_now;
				}
				
				// Now write out audio.
				audioWrite(context, n, channel, out);
			}
			
		//gScope.log(out_freqs[0], out_freqs[1]);
		
	}

}

void cleanup(BelaContext *context, void *userData)
{

}
