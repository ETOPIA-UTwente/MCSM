import time, os, pickle, threading, argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import butter, lfilter, freqz, blackmanharris
import math
#from picosdk.functions import adc2mV


status = {}
parser = argparse.ArgumentParser()
parser.add_argument('--device_list', required=True)
parser.add_argument('--start_time_offset', required=False, default=0)
parser.add_argument('--test_name', required=True)

args = parser.parse_args()
device_list = args.device_list.split(",")                                                     # in seconds [ Total time to stream]
start_time_offset = float(args.start_time_offset)
test_name = args.test_name
pol_shift = True
start_time_offset = 0

# Where the data is stored for a device
file_names = {}
data_directories = {}
sample_period = 0
for i, device in enumerate(device_list):
  if device != "":
    data_directories[device] = f"C:\\pi_data\\{device}\\{test_name}\\dataset"
    file_names[device] = data_directories[device]

print(data_directories)
print(device_list)

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def filter(data, n, order, fs, cutoff):
    # Filter requirements.
    """ order = 6
    fs = 30.0       # sample rate, Hz
    cutoff = cutoff_freq  # desired cutoff frequency of the filter, Hz """

    # Get the filter coefficients so we can check its frequency response.
    b, a = butter_lowpass(cutoff, fs, order)
    # Plot the frequency response.
    w, h = freqz(b, a, worN=8000)
    t = np.linspace(0, n , n, endpoint=False)
    # Filter the data, and plot both the original and filtered signals.
    y = butter_lowpass_filter(data, cutoff, fs, order)
    return y

data = {}
params = []
temp_data = []



start_time_list = {}
key_map_scale = {}
key_map_scale["pi_1"] = { "CH1":  100,           # V Res
                          "CH2":  100,           # V com 
                          "CH3":  1000,          # I Leds
                          "CH4":  100,          # I Water pump
                          "CH5":  1000,           # I Heater
                          "CH6":  1000,           # I Drill
                          "trigger": 1000,
                          "trigger_index": 1 
                        }
key_map_scale["pi_2"] = { "CH1": 100,           # V supply
                          "CH2": 100,           # I Supply
                          "CH3": 1000,           # I VFD
                          "CH4": 1000,
                          "CH5": 100,
                          "CH6": 100,
                          "trigger": 1000,
                          "trigger_index": 1 
                        }


for device in data_directories:
  data[device] = {  "CH1": np.array([]),
                    "CH2": np.array([]),
                    "CH3": np.array([]),
                    "CH4": np.array([]),
                    "CH5": np.array([]),
                    "CH6": np.array([]),
                    "trigger": np.array([]),
                    "trigger_index": np.array([])
                    }
          
  read = open(f"{file_names[device]}", "rb")
  temp_data = pickle.load(read, encoding="bytes")

  """ print("Max ----------------------- \n")
  print(max(temp_data["CH1"]))
  print(temp_data["max_adc"])
  print("----------\n") """

  read.close()

  v_scale = 100
  i_scale = 100
  ch_scale = 1000
  
  sample_period = float(temp_data["speriod"])
  sample_freq = 1 / sample_period
  print(f"Sample Freq : {sample_freq}")
  """ class adc:
    def __init__(self, value):
      self.value = value
    
  adc = adc(temp_data["max_adc"])
  adc2mVChAMax = adc2mV(temp_data["CH1"], 8, adc)

  print("Max ----------------------- \n")
  print(max(adc2mVChAMax))
  #print(temp_data["max_adc"])
  print("----------\n") """

  #adc_multiplier = (5000 / int(temp_data["max_adc"]))

  #how many data points
  number_data_points = len(temp_data["CH1"]) #  500e6
  #scale data from transducers - 
  for key in data[device]:
    max_adc_scaled = temp_data["max_adc"] / 10#key_map_scale[device][key]
    adc_multiplier = (5000 / max_adc_scaled)
    if len(temp_data[key]) <= 500e6:  
      data[device][key] = np.append(data[device][key], filter(temp_data[key], len(temp_data[key]), 5, int(sample_freq), int(sample_freq)*0.4))
      data[device][key] = (data[device][key] / key_map_scale[device][key]) * adc_multiplier
    else:
      data[device][key] = np.append(data[device][key], filter(temp_data[key][0:number_data_points], len(temp_data[key]), 5, int(sample_freq), int(sample_freq)*0.4))
      data[device][key] = (data[device][key] / key_map_scale[device][key]) * adc_multiplier

  data[device]["trigger_time"] = float(temp_data["trigger_time"])




pi_1_trig_time = float(data["pi_1"]["trigger_time"])
pi_2_trig_time = float(data["pi_2"]["trigger_time"]) - float(start_time_offset)
if "pi_3" in device_list:
  pi_3_trig_time = float(data["pi_3"]["trigger_time"]) #- float(start_time_offset)

pi_1_trig_index = int(data["pi_1"]["trigger_index"])
pi_2_trig_index = int(data["pi_2"]["trigger_index"])
print(f"Pi_1 trig index : {pi_1_trig_index} Pi_2 trig index: {pi_2_trig_index}")
if "pi_3" in device_list:
  pi_3_trig_index = int(data["pi_3"]["trigger_index"])

print(f"pi1 trig time - {pi_1_trig_time}")
print(f"pi2 trig time - {pi_2_trig_time}")
if "pi_3" in device_list:
  print(f"pi3 trig time - {pi_3_trig_time}")
print(f"Sample period - {sample_period}")

min_index = min(pi_1_trig_index, pi_2_trig_index)#, pi_3_trig_index)
sync_pi_1 = pi_1_trig_index - min_index #sorted_index[sorted_index.keys().index("pi_1")
sync_pi_2 = pi_2_trig_index - min_index
if "pi_3" in device_list:
  sync_pi_3 = pi_3_trig_index - min_index

shift_1_2 = 0
if "pi_3" in device_list:
  shift_1_3 = 0
  pol_shift_1_3 = 0
pol_shift_1_2 = 0

if pol_shift == True:
  sync_data_len = len(data["pi_1"]["CH1"]) 
  if len(data["pi_1"]["CH1"]) > (3 / sample_period):
    sync_data_len = int(3 / sample_period)
    print("reducing length of sync data to 3 / sample_period")

  sync_V1 = data["pi_1"]["CH1"][sync_pi_1:sync_data_len + sync_pi_1][0:sync_data_len]
  print("\n\n --------------")
  print(len(sync_V1))
  sync_V1 = -filter(sync_V1, len(sync_V1), 5, 1250000, 1000)
  sync_V2 = data["pi_2"]["CH1"][sync_pi_2:sync_data_len + sync_pi_2][0:sync_data_len]
  sync_V2 = -filter(sync_V2, len(sync_V2), 5, 1250000, 1000)
  print(len(sync_V2))
  if "pi_3" in device_list:
    sync_V3 = data["pi_3"]["CH1"][sync_pi_3: sync_data_len + sync_pi_3][0:sync_data_len]
    sync_V3 = -filter(sync_V3, len(sync_V3), 5, 1250000, 1000)

  #data shift algo   
  try:
    #a = 1 / 0
    phi_1_2 = math.acos( np.dot(sync_V1,sync_V2) / (np.linalg.norm(sync_V1)*np.linalg.norm(sync_V2)) )
    phi_1_2_deg = (phi_1_2 / math.pi) * 180

    if "pi_3" in device_list:
      print(f"Phi 1-2: {phi_1_2} - {phi_1_2_deg}deg")
      phi_1_3 = math.acos( np.dot(sync_V1,sync_V3) / (np.linalg.norm(sync_V1)*np.linalg.norm(sync_V3)) )
      phi_1_3_deg = (phi_1_3 / math.pi) * 180
      print(f"Phi 1-3: {phi_1_3} - {phi_1_3_deg}deg")

    one_rev = 20e-3 / sample_period
    print(f"One rev = {one_rev}")
    shift_1_2 = int((phi_1_2_deg / 360) * one_rev)
    if "pi_3" in device_list:
      shift_1_3 = int((phi_1_3_deg / 360) * one_rev)
      max_shift = max(abs(shift_1_2), abs(shift_1_3))
    else:
      max_shift = abs(shift_1_2)

    pol_sync = int(sync_data_len*0.5)
    pol_sync_V1 = sync_V1[max_shift: max_shift + pol_sync]
    pol_sync_V2 = sync_V2[max_shift - shift_1_2 : pol_sync + max_shift - shift_1_2 ]
    if "pi_3" in device_list:
      pol_sync_V3 = sync_V2[max_shift - shift_1_3 : pol_sync + max_shift - shift_1_3 ]

    pol_phi_1_2 = math.acos( np.dot(pol_sync_V1,pol_sync_V2) / (np.linalg.norm(pol_sync_V1)*np.linalg.norm(pol_sync_V2)) )
    pol_phi_1_2_deg = int((pol_phi_1_2 / math.pi) * 180)

    if "pi_3" in device_list:
      pol_phi_1_3 = math.acos( np.dot(pol_sync_V1,pol_sync_V3) / (np.linalg.norm(pol_sync_V1)*np.linalg.norm(pol_sync_V3)) )
      pol_phi_1_3_deg = int((pol_phi_1_3 / math.pi) * 180)
      pol_shift_1_3 = int((pol_phi_1_3_deg / 360) * one_rev)

    pol_shift_1_2 = int((pol_phi_1_2_deg / 360) * one_rev)
  except Exception as err:
    print(f"Could not do phase correction - Output may be DC Error: {err}")
    

  polarity_1_2 = 1
  if "pi_3" in device_list:
    polarity_1_3 = 1
  if pol_shift_1_2 > shift_1_2:
    polarity_1_2 = -1

  if "pi_3" in device_list:
    if pol_shift_1_3 > shift_1_3:
      polarity_1_3 = -1

  print(f"shift 1 - 2 = {polarity_1_2*shift_1_2}")
  if "pi_3" in device_list:
    print(f"shift 1 - 3 = {polarity_1_3*shift_1_3}")
  print(f"POL shift 1 - 2 = {pol_shift_1_2}")
  if "pi_3" in device_list:
    print(f"POL shift 1 - 3 = {pol_shift_1_3}")

  max_shift = abs(shift_1_2) # max(abs(shift_1_2), abs(shift_1_3))

  #max_shift = 0
  sync_pi_1 = sync_pi_1 + max_shift 
  sync_pi_2 = sync_pi_2 + max_shift - shift_1_2*polarity_1_2
  if "pi_3" in device_list:
    sync_pi_3 = sync_pi_3 + max_shift - shift_1_3*polarity_1_3



print(f"pi1 trig index - {pi_1_trig_index}")
print(f"pi2 trig index - {pi_2_trig_index}")
if "pi_3" in device_list:
  print(f"pi3 trig index - {pi_3_trig_index}")
print("-------------------------------")
print(f"sync pi 1 index - {sync_pi_1}")
print(f"sync pi 2 index - {sync_pi_2}")
if "pi_3" in device_list:
  print(f"sync pi 3 index - {sync_pi_3}")



data_len_pi_1 = len(data["pi_1"]["CH1"])
data_len_pi_2 = len(data["pi_2"]["CH1"]) 
if "pi_3" in device_list:
  data_len_pi_3 = len(data["pi_3"]["CH1"]) 

if "pi_3" in device_list:
  time_sync_len = min(data_len_pi_1, data_len_pi_2, data_len_pi_3) - max(sync_pi_1, sync_pi_2, sync_pi_3) 
else:
  time_sync_len = min(data_len_pi_1, data_len_pi_2) - max(sync_pi_1, sync_pi_2)

print("Loading complete")
start_index = 0
end_index_ref = time_sync_len  #- 1000000#int(len(data["pi_1"]["voltage_L1"]))# 500000

print(end_index_ref)
while start_index < end_index_ref:
  end_index = start_index + end_index_ref - 1000000
  if end_index >= end_index_ref:
    end_index = end_index_ref

  end_index = end_index_ref


  time = np.linspace(start_index, end_index, end_index - start_index)
  time = time * sample_period * 1000 # set x-axis to milliseconds
  print("\n\nsample period")
  print(sample_period)
  print(data["pi_1"].keys())
  print(data["pi_2"].keys())
  if "pi_3" in device_list:
    print(data["pi_3"].keys())

  

  

  

  V1 = data["pi_1"]["CH1"][start_index + sync_pi_1: end_index + sync_pi_1]
  V1 = -filter(V1, len(V1), 5, 1250000, 125000)
  V2 = data["pi_2"]["CH1"][start_index + sync_pi_2: end_index + sync_pi_2]
  V2 = -filter(V2, len(V2), 5, 1250000, 125000)
  if "pi_3" in device_list:
    V3 = data["pi_3"]["CH1"][start_index + sync_pi_3: end_index + sync_pi_3]
    V3 = -filter(V3, len(V3), 10, 1250000, 20000)

  I1 = data["pi_1"]["CH2"][start_index + sync_pi_1: end_index + sync_pi_1]
  I1 = filter(I1, len(I1), 5, 1250000, 125000)
  I2 = data["pi_1"]["CH2"][start_index + sync_pi_1: end_index + sync_pi_1]
  I2 = filter(I2, len(I2), 5, 1250000, 125000)
  I3 = data["pi_2"]["CH2"][start_index + sync_pi_2: end_index + sync_pi_2]
  I3 = filter(I3, len(I3), 7, 1250000, 10000)


  print("Offset current 1,2,3")
  print( np.mean(I1[0:1000000], axis=0))
  I1_mean_offset = np.mean(I1[0:1000000], axis=0)
  print( np.mean(I2[0:1000000], axis=0))
  I2_mean_offset = np.mean(I2[0:1000000], axis=0)
  print( np.mean(I3[0:1000000], axis=0))
  I3_mean_offset = np.mean(I3[0:1000000], axis=0)

  

  

  """ try:
    begin_index, end_index = 3800000,7600000
    begin_index, end_index = 2000000,6600000  
    #RMS Current calculation
    ms = 0
    N = 0
    for i in range(begin_index, end_index):
        ms = ms + I1[i]*I1[i]
        N += 1
    ms = ms / N
    I_rms_1 = math.sqrt(ms) - I1_mean_offset
    print(f"RMS current (1): {I_rms_1}")

    ms = 0
    N = 0
    for i in range(begin_index, end_index):
        ms = ms + I2[i]*I2[i]
        N += 1
    ms = ms / N
    I_rms_2 = math.sqrt(ms) - I2_mean_offset
    print(f"RMS current (2): {I_rms_2}")
    ms = 0
    N = 0
    for i in range(begin_index, end_index):
        ms = ms + I3[i]*I3[i]
        N += 1
    ms = ms / N
    I_rms_3 = math.sqrt(ms) - I3_mean_offset
    print(f"RMS current (3): {I_rms_3}")
    #RMS Current calculation
    ms = 0
    N = 0
    for i in range(begin_index, end_index):
        ms = ms + V1[i]*V1[i]
        N += 1
    ms = ms / N
    V_rms_1 = math.sqrt(ms)
    print(f"RMS Voltage (1): {V_rms_1}")

    ms = 0
    N = 0
    for i in range(begin_index, end_index):
        ms = ms + V2[i]*V2[i]
        N += 1
    ms = ms / N
    V_rms_2 = math.sqrt(ms)
    print(f"RMS voltage (2): {V_rms_2}")
    ms = 0
    N = 0
    for i in range(begin_index, end_index):
        ms = ms + V3[i]*V3[i]
        N += 1
    ms = ms / N
    V_rms_3 = math.sqrt(ms) #- I3_mean_offset
    print(f"RMS Voltage (3): {V_rms_3}")

    print(f"RMS Power (1): {V_rms_1*I_rms_1}")
    print(f"RMS Power (2): {V_rms_2*I_rms_2}")
    print(f"RMS Power (3): {V_rms_3*I_rms_3}")
  except Exception as err:
    print(f"Exception occured {err}") """

  


  #RMS voltage
  """ V = V1
  I = I1
  #rms = norm(V)/sqrt(V.size)
  V_rms = np.sqrt(np.mean(V**2))
  I_rms = np.sqrt(np.mean(I**2)) """
  #rms = lambda V, axis=None: np.sqrt(np.mean(np.square(V), axis))
  """ print(f"RMS voltage: {V_rms}")
  print(f"RMS current: {I_rms}")
  print(f"RMS power: {I_rms*V_rms}") """

  """ power_1 = V1 * I1
  power_2 = V2 * I2
  power_3 = V3 * I3 """


  #offset current correction I3
  #I3 = I3 - I3_mean_offset


  start_point = 0
  end_point = len(data["pi_1"]["CH1"])#500000

  
  y_temp = np.linspace(0, end_point - start_point, end_point - start_point)
  y_temp = y_temp * sample_period * 1000 # set x-axis to milliseconds
  
  plt.figure(1)
  plt.plot(y_temp, data["pi_1"]["CH3"][start_point: end_point], color="blue", label="I_(LEDs)", linewidth=2)
  plt.plot(y_temp, data["pi_1"]["CH4"][start_point: end_point], color="red", label="I_(pump)", linewidth=2)
  plt.plot(y_temp, data["pi_1"]["CH5"][start_point: end_point], color="green", label="I_(heater)", linewidth=2)
  plt.plot(y_temp, data["pi_1"]["CH6"][start_point: end_point], color="orange", label="I_(drill)", linewidth=2)

  plt.legend(loc='upper right')
  plt.hlines(y=0.0, xmin=0, xmax=40, linewidth=1, color='black')
  plt.ylabel('Current (A)')
  plt.xlabel('Time (ms)')
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)
  plt.gca().spines['left'].set_linewidth(2)
  plt.gca().spines['bottom'].set_linewidth(1.5)
  plt.ylim((-100, 100))
  
  plt.figure(2)
  plt.plot(y_temp, data["pi_2"]["CH1"][start_point: end_point], color="blue", label="V(Main)", linewidth=2)
  plt.plot(y_temp, data["pi_1"]["CH1"][start_point: end_point], color="red", label="V(Res)", linewidth=2)
  plt.plot(y_temp, data["pi_1"]["CH2"][start_point: end_point], color="green", label="V(Com)", linewidth=2)

  plt.legend(loc='upper right')
  plt.hlines(y=0.0, xmin=0, xmax=40, linewidth=1, color='black')
  plt.ylabel('Voltage (V)')
  plt.xlabel('Time (ms)')
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)
  plt.gca().spines['left'].set_linewidth(2)
  plt.gca().spines['bottom'].set_linewidth(1.5)
  plt.ylim((-400, 400))

  plt.figure(3)  
  plt.plot(y_temp, data["pi_2"]["CH2"][start_point: end_point], color="blue", label="I_supply", linewidth=2)
  #plt.plot(y_temp, -data["pi_2"]["CH3"][start_point: end_point], color="green", label="I(Supply - ch3)", linewidth=2)
  #plt.plot(y_temp, -data["pi_2"]["CH4"][start_point: end_point], color="red", label="I(Supply - ch4)", linewidth=2)
  #plt.plot(y_temp, data["pi_2"]["CH3"][start_point: end_point], color="red", label="I(VFD)", linewidth=2)
 
  plt.legend(loc='upper right')
  plt.hlines(y=0.0, xmin=0, xmax=40, linewidth=1, color='black')
  plt.ylabel('Current (A)')
  plt.xlabel('Time (ms)')
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)
  plt.gca().spines['left'].set_linewidth(2)
  plt.gca().spines['bottom'].set_linewidth(1.5)
  plt.ylim((-200, 200))
  

  plt.figure(4)
  plt.plot(time, V1, color="blue", label="V_(Res)", linewidth=2)
  plt.plot(time, V2, color="red", label="V_(Supply)", linewidth=2)

  plt.legend(loc='upper right')
  plt.hlines(y=0.0, xmin=0, xmax=40, linewidth=1, color='black')
  plt.ylabel('Voltage (V)')
  plt.xlabel('Time (ms)')
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)
  plt.gca().spines['left'].set_linewidth(2)
  plt.gca().spines['bottom'].set_linewidth(1.5)
  plt.ylim((-400, 400))

  

  #instantaneous power
  plt.figure(6)
  #total_current = data["pi_1"]["CH3"][start_point: end_point] + data["pi_1"]["CH4"][start_point: end_point] + data["pi_1"]["CH5"][start_point: end_point] + data["pi_1"]["CH6"][start_point: end_point]
  supply_current = data["pi_2"]["CH2"]
  supply_voltage = data["pi_2"]["CH1"]
  #plt.plot(y_temp, abs(total_current * data["pi_1"]["CH1"][start_point: end_point])  , color="blue", label="P_(total)", linewidth=2)
  plt.plot(y_temp, abs(supply_current * supply_voltage)  , color="blue", label="P_(total)", linewidth=2)
  #plt.plot(y_temp, supply_voltage  , color="red", label="V_(total)", linewidth=2)
  #plt.plot(y_temp, supply_current, color="green", label="I_(total)", linewidth=2)
  #plt.plot(y_temp, data["pi_1"]["CH4"][start_point: end_point], color="red", label="I_(drill)", linewidth=2)
  #plt.plot(y_temp, data["pi_1"]["CH5"][start_point: end_point], color="green", label="I_(heater)", linewidth=2)
  #plt.plot(y_temp, data["pi_1"]["CH6"][start_point: end_point], color="orange", label="I_(pump)", linewidth=2)
  plt.legend(loc='upper right')
  plt.hlines(y=0.0, xmin=0, xmax=40, linewidth=1, color='black')
  plt.ylabel('Instantaneous Power (W)')
  plt.xlabel('Time (ms)')
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)
  plt.gca().spines['left'].set_linewidth(2)
  plt.gca().spines['bottom'].set_linewidth(1.5)
  plt.ylim((-1000, 30000))

  power_data_set = abs(supply_current * supply_voltage)
  inrush_index_tracker = []
  inrush_zero_tracker = []
  data_incomplete = True
  index = 0
  while (data_incomplete):
    try:
      x = power_data_set[index]
    except:
      data_incomplete = False
      break
    if power_data_set[index] > 5500:
      inrush_index_tracker.append(index)
      inrush_zero_tracker.append(0)
      index += 180000
    index += 1
  string_index = "["
  for item in inrush_index_tracker:
    string_index += f"{item}], ["
  string_index += "]"
  print(string_index)
  print(inrush_zero_tracker)
  

  #instantaneous power
  """ plt.figure(7)
  total_current = data["pi_1"]["CH3"][start_point: end_point] + data["pi_1"]["CH4"][start_point: end_point] + data["pi_1"]["CH5"][start_point: end_point] + data["pi_1"]["CH6"][start_point: end_point]
  plt.plot(y_temp, abs(total_current * data["pi_1"]["CH1"][start_point: end_point])  , color="blue", label="P_(total)", linewidth=2)

  plt.legend(loc='upper right')
  plt.hlines(y=0.0, xmin=0, xmax=40, linewidth=1, color='black')
  plt.ylabel('Instantaneous Power (W)')
  plt.xlabel('Time (ms)')
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)
  plt.gca().spines['left'].set_linewidth(2)
  plt.gca().spines['bottom'].set_linewidth(1.5)
  plt.ylim((-10000, 10000)) """

  fig, ax1 = plt.subplots()
  color = 'tab:red'
  ax1.set_xlabel('Time (ms)', fontsize=10)
  ax1.set_ylabel('Current (A)', color=color, fontsize=10)
  ax1.plot(y_temp, supply_current, color=color, linewidth=2)
  ax1.tick_params(axis='y', labelcolor=color, labelsize=10)
  ax1.set_ylim((-100, 100))

  ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

  color = 'tab:blue'
  ax2.set_ylabel('Voltage (V)', color=color, fontsize=10)  # we already handled the x-label with ax1
  ax2.plot(y_temp, -supply_voltage, color=color, linewidth=2)
  ax2.tick_params(axis='y', labelcolor=color, labelsize=10)
  ax2.set_ylim((-400, 400))

  """ plt.legend(loc='upper right')
  plt.hlines(y=0.0, xmin=0, xmax=40, linewidth=1, color='black')
  plt.ylabel('Power (W)')
  plt.xlabel('Time (ms)')
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)
  plt.gca().spines['left'].set_linewidth(2)
  plt.gca().spines['bottom'].set_linewidth(1.5)
  plt.ylim((-1000, 20000)) """
  


  plt.show()
  start_index = end_index