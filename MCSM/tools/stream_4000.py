# Copyright (C) 2018-2019 Pico Technology Ltd. See LICENSE file for terms.
# PS2000 Series (A API) STREAMING MODE EXAMPLE
# This example demonstrates how to call the ps4000 driver API functions in order to open a device, setup 2 channels and collects streamed data (1 buffer).
# This data is then plotted as mV against time in ns.
#########################[ Stream for amount of time [ seconds ] ]###############################(Amatthee-22-10-2019)####################
# Description: I have modified this script to run with a set time window in seconds. The original script ran a set specified buffersize  #
##########################################################################################################################################
import PI_Functions 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--time', required=True)
parser.add_argument('--speriod', required=True)
parser.add_argument('--tbc', required=False)
parser.add_argument('--timeStamp', required=False)
parser.add_argument('--rpi', required=False)
parser.add_argument('--start_time', required=False)
args = parser.parse_args()
# First step is to sync th epi clock
if args.timeStamp is not None and args.rpi is not None:
  raspberry = PI_Functions.Setup_PI(float(args.timeStamp))
  raspberry.setup_GPIO(21)

import ctypes, time, threading, traceback, argparse, psutil, pickle, os, sys, math
from picosdk.functions import adc2mV, assert_pico_ok
from admin_4000 import Admin
import numpy as np

try:
  # Create chandle and status ready for use
  status = {}
  streamTimeWindow = float(args.time)                                                      # in seconds [ Total time to stream]
  sampleTimeInNS = int(args.speriod) 
  try:
    timeBetweenCallbacks = float(args.tbc)
  except:
    timeBetweenCallbacks = None
  try:
    ref_time = float(args.timeStamp)
  except:
    ref_time = time.time()
  try:
    start_time = float(args.start_time)
  except:
    start_time = ref_time + 10
  time_offset = 0
  rpi = True
  #########################[ Stream for amount of time [ seconds ] ]###############################(Amatthee-22-10-2019)####################
  # Description: I have modified this script to run with a set time window in seconds. The original script ran a set specified buffersize  #
  ##########################################################################################################################################
  # Size of capture # previous buffer size code
  realSampleTime = sampleTimeInNS*1e-9
  nextSample = 0
  enumerator = 0
  startIndex, sourceEnd = 0,0
  trigger_time = 0
  number_of_channels = 4
  process = Admin(number_of_channels, realSampleTime, raspi=rpi)
  append_file = False
  process.sizeOfOneBuffer = int(streamTimeWindow / realSampleTime)
  if timeBetweenCallbacks:
    process.sizeOfOneBuffer = int(timeBetweenCallbacks / realSampleTime)

  callbacks_buff_len = process.sizeOfOneBuffer
  callbacks_buff_len_counter = 0
  totalSamples = int(streamTimeWindow / realSampleTime)
  print(f"total samples - {totalSamples}")
  print("#"*80)
  print(f"#\t Start time (time stamp) : {ref_time}")
  print(f"#\t Max buffer to read (samples): {process.sizeOfOneBuffer}")
  print(f"#\t Streaming time period (S): {streamTimeWindow}")
  print(f"#\t Sample period (nS): {sampleTimeInNS} \t sample freq: {1/realSampleTime} HZ")
  print(f"#\t Time between callbacks (S): {timeBetweenCallbacks}")
  print(f"#\t Append file function is : {append_file}")
  print("#"*80)
  ##########################################################################################################################################
  process.picoSetupHandle()
  channels = [0,1,2,3,4,5,6]
  channels_ = ["A", "B", "C", "D", "E", "F", "G"]
  for channel in channels:
      process.picoSetupChannel(channel)

  process.time_units(sampleTimeInNS)
  process.totalSamples = totalSamples
  trigger_complete = False

  bufferAMax = np.zeros(shape=process.sizeOfOneBuffer, dtype=np.int16)
  bufferBMax = np.zeros(shape=process.sizeOfOneBuffer, dtype=np.int16)
  bufferCMax = np.zeros(shape=process.sizeOfOneBuffer, dtype=np.int16)
  bufferDMax = np.zeros(shape=process.sizeOfOneBuffer, dtype=np.int16)
  bufferEMax = np.zeros(shape=process.sizeOfOneBuffer, dtype=np.int16)
  bufferFMax = np.zeros(shape=process.sizeOfOneBuffer, dtype=np.int16)
  bufferGMax = np.zeros(shape=process.sizeOfOneBuffer, dtype=np.int16)

  #for channel in channels:
  process.picoSetupChannelBuffer(channels_[0], bufferAMax)
  process.picoSetupChannelBuffer(channels_[1], bufferBMax)
  process.picoSetupChannelBuffer(channels_[2], bufferCMax)
  process.picoSetupChannelBuffer(channels_[3], bufferDMax)
  process.picoSetupChannelBuffer(channels_[4], bufferEMax)
  process.picoSetupChannelBuffer(channels_[5], bufferFMax)
  process.picoSetupChannelBuffer(channels_[6], bufferGMax)
  
  waitTimeSec = 10
  # Convert the python function into a C function pointer.
  condition = True # condition to continue streaming
  start_object = math.trunc(ref_time)

  start_condition = True
  print(f"waiting for {waitTimeSec} sec - starting at {start_time}")
  print(f"time now is {time.time()}")
  
  while start_condition:
    #process.ps.ps4000aFlashLed(process.chandle, 1)
    #start_condition = (time.time() < start_object + waitTimeSec - time_offset)  # wait for 2 sec
    start_condition = (time.time() < start_time)  # wait for 2 sec
  
  #Trigger has a delay to allow for deviations in run_stream call
  trigger_time = time.time()
  process.setTrigger(0)
  delay = 1
  input_port = 21
  no_pulse_edges = 4
  pulse_width = 0.5
  pulse_width_decay = 0.1
  trig = threading.Thread(name='child procs', target=process.fireTrigger, args=(input_port, delay, no_pulse_edges, pulse_width, pulse_width_decay, trigger_time))
  trig.start()

  process.run_streaming()
  
  # Caculate some useless info...
  actualSampleInterval = process.sampleInterval.value
  actualSampleIntervalNs = actualSampleInterval 
  actualSampleIntervalScaled = actualSampleIntervalNs
  actualSampleIntervalScaledUnit = "ns" 
  # Print out sample interval info in correct units
  if actualSampleIntervalNs >= 1000:
      actualSampleIntervalScaled = actualSampleIntervalNs / 1000
      actualSampleIntervalScaledUnit = "us"
      if actualSampleIntervalScaled >= 1000:
          actualSampleIntervalScaled = actualSampleIntervalScaled / 1000
          actualSampleIntervalScaledUnit = "ms"

  print("#"*80)
  print(f"# \tCapturing at sample interval {actualSampleIntervalScaled} {actualSampleIntervalScaledUnit}")
  print(f"# \tSampling window time is {streamTimeWindow} S. Total sample dataset of {totalSamples} samples.")
  print(f"# \tTime between callbacks {timeBetweenCallbacks}")
  print("#"*80)

  def enumerateBinaryFile(enumerator, bufferAMax, bufferBMax, bufferCMax, bufferDMax, bufferEMax, bufferFMax, bufferGMax):
    global start_time, realSampleTime, trigger_time
    file_name = f"/home/pi/Desktop/pi_measure/code/data/"
    while int(psutil.virtual_memory()._asdict()['percent']) > 70:
      print("Thread: Im waiting for some memory...")
      time.sleep(5) # wait for processes
    data = {
          "start_time"    :  start_time,
          "trigger_time"  :  trigger_time,
          "speriod"       :  realSampleTime,
          "max_adc"       :  process.maxADC,
          "CH1"            :  bufferAMax,
          "CH2"            :  bufferBMax,
          "CH3"            :  bufferCMax,
          "CH4"            :  bufferEMax,
          "CH5"            :  bufferFMax,
          "CH6"            :  bufferGMax,
          "trigger"        :  bufferDMax#,
          #"time"          :  
          }
    pickle_out = open(f"{file_name}{enumerator}","wb")
    pickle.dump(data, pickle_out)
    pickle_out.close()
    data = None

  def add_trigger_index(trigger_index):
    file_name = f"/home/pi/Desktop/pi_measure/code/data/"
    while int(psutil.virtual_memory()._asdict()['percent']) > 70:
      print("Thread: Im waiting for some memory...")
      time.sleep(5) # wait for processes
    
    file0 = open(f"{file_name}0", "rb")
    master_data = {}

    file_0 = pickle.load(file0)

    for key in file_0.keys():
      master_data[key] = np.array([])
      master_data[key] = np.append(master_data[key], file_0[key])
    master_data["trigger_index"] = np.array([trigger_index])
    file0.close()
    pickle_out = open(f"{file_name}0","wb")
    pickle.dump(master_data, pickle_out)
    pickle_out.close()

  def streaming_callback(handle, noOfSamples, startIndex, overflow, triggerAt, triggered, autoStop, param):
      global nextSample, autoStopOuter, wasCalledBack, sampleStartTime, enumerator, trigger_index
      global trigger_flag, callbacks_buff_len_counter, callbacks_buff_len 
      global buff_to_store_a, buff_to_store_b, buff_to_store_c, buff_to_store_d, buff_to_store_e, buff_to_store_f, buff_to_store_g

      wasCalledBack = True
      if triggered == 1 or triggered == '1':
        print(f"Triggered at ----- {triggerAt}")
        trigger_index += triggerAt
        trigger_flag = False
      if trigger_flag:
        trigger_index += noOfSamples
      sampleStartTime = time.time()
      destEnd = nextSample + noOfSamples
      sourceEnd = startIndex + noOfSamples
      #barCount += noOfSamples
      nextSample += noOfSamples
      callbacks_buff_len_counter = callbacks_buff_len_counter + noOfSamples

      buff_to_store_a = np.append(buff_to_store_a, bufferAMax[startIndex:sourceEnd])
      buff_to_store_b = np.append(buff_to_store_b, bufferBMax[startIndex:sourceEnd])
      buff_to_store_c = np.append(buff_to_store_c, bufferCMax[startIndex:sourceEnd])
      buff_to_store_d = np.append(buff_to_store_d, bufferDMax[startIndex:sourceEnd])
      buff_to_store_e = np.append(buff_to_store_e, bufferEMax[startIndex:sourceEnd])
      buff_to_store_f = np.append(buff_to_store_f, bufferFMax[startIndex:sourceEnd])
      buff_to_store_g = np.append(buff_to_store_g, bufferGMax[startIndex:sourceEnd])
     
      if callbacks_buff_len_counter >= callbacks_buff_len:
        t = threading.Thread(name='child procs', target=enumerateBinaryFile, args=(enumerator,buff_to_store_a, 
                                                                                              buff_to_store_b, 
                                                                                              buff_to_store_c, 
                                                                                              buff_to_store_d,
                                                                                              buff_to_store_e,
                                                                                              buff_to_store_f,
                                                                                              buff_to_store_g))
        t.start()
        callbacks_buff_len_counter = 0
        enumerator += 1
        buff_to_store_a = np.array([])
        buff_to_store_b = np.array([])
        buff_to_store_c = np.array([])
        buff_to_store_d = np.array([]) 
        buff_to_store_e = np.array([])
        buff_to_store_f = np.array([])
        buff_to_store_g = np.array([])
      
  buff_to_store_a = np.array([])
  buff_to_store_b = np.array([])
  buff_to_store_c = np.array([])
  buff_to_store_d = np.array([]) 
  buff_to_store_e = np.array([])
  buff_to_store_f = np.array([])
  buff_to_store_g = np.array([])
  # Convert the python function into a C function pointer.
  cFuncPtr = process.ps.StreamingReadyType(streaming_callback)
  # Fetch data from the driver in a loop, copying it out of the registered buffers and into our complete one.
  condition = True # condition to continue streaming
  autoStopOuter = False
  nextSample = 0
  enumerator = 0
  trigger_index = 0
  trigger_flag = True
  start_time = time.time() # start time for streaming
  sampleStartTime = start_time

  print(f"Sample start time {sampleStartTime}")
  while condition:
      wasCalledBack = False
      process.get_latests_stream_values(cFuncPtr)
      condition = nextSample < totalSamples and not autoStopOuter
      if not wasCalledBack:
        time.sleep(0.01)

  print(f"Trigger Index -*-{trigger_index}-*-")
  add_trigger_index(trigger_index)
  end_time = time.time()
  total_time = end_time - start_time
  process.stop_and_close()
except Exception as err:
  print(f"Error: {err}")
  traceback.print_exc()
  process.stop_and_close()

