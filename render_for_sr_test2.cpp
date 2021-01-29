```
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
	
	std::vector<float> logs(1000000);
	unsigned int writeFileSize = logs.size() + 5000;
	freq_mod_tone_bela.setBufferSize(writeFileSize);
	unsigned int logIdx = 0;
	
	while(!Bela_stopRequested())
	{
		
		float data_out_now[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		 data_out_now[0] = { timestamp_ms };
		
		for(unsigned int n1 = 0; n1 < 10; ++n1)
				{
					logs[logIdx++] = data_out_now[n1]; // fill an array in memory

						if(logs.size() == logIdx)
						{
							// when the array is full, dump it to disk
							// given how we have setBufferSize() above, WriteFile's buffer
							// should have enough space for all of logs
							printf("Dumping\n");
							freq_mod_tone_bela.log(logs.data(), logIdx);
							logIdx = 0;
							// now we need to wait until the WriteFile thread has written everything to disk
							// unfortunately it will always keep up to 4096 samples inside its buffer
							while(freq_mod_tone_bela.getBufferStatus() < (writeFileSize - logs.size()) / writeFileSize)
							{
								printf("waiting %f\n", freq_mod_tone_bela.getBufferStatus());
								usleep(10000);
							}
							printf("Running\n");
						}
					}

		// Read locations from Trill sensor
		unsigned int q = touchSensor.size() < NUM_TOUCH ? touchSensor.size() : NUM_TOUCH;
		for(unsigned int n = 0; n < q; ++n){

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
				if(read_out[n]!=past_read[n]){
					//float data_out_now[5] =	{ timestamp_ms, float(n), read_out[n], y_coordinate_now[n], frequency[n] };
					//instead of logging every time ...
					//freq_mod_tone_bela.log(data_out_now,5);
					//we log a larger buffer at once, so the WriteFile thread does not run in parallel to us
					
				}
			} else {
				active_touch[n]=0.0;
			}
		}

	//usleep(gTaskSleepTime);
	}
	// send WriteFile any leftover stuff
	freq_mod_tone_bela.log(logs.data(), logIdx);
}


bool setup(BelaContext *context, void *userData)
{
	
	
	// output file.
	freq_mod_tone_bela.setup("freq_mod_tone_bela.bin");
	freq_mod_tone_bela.setFormat("%f %f %f %f %f %f %f %f %f %f\n");
	freq_mod_tone_bela.setFileType(kBinary); 
	//freq_mod_tone_bela.setBufferSize(writeFileSize);
	
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

```
