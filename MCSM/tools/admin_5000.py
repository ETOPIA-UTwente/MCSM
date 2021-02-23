import os, ctypes, pprint, pylab as pl, time, datetime
#from tools import database
from picosdk.ps5000a import ps5000a as ps
from picosdk.functions import adc2mV, mV2adc, assert_pico_ok
from ctypes import Structure
import RPi.GPIO as GPIO
#import monotonic
#mtime = monotonic.time.time
#from picosdk.ps5000a import PS5000A_TRIGGER_INFO

pp = pprint.PrettyPrinter(indent=4, width=80)

class Admin(object):
  class PS5000A_TRIGGER_INFO(Structure):
    _pack_ = 1
    _fields_ = [("status", ctypes.c_uint32),
                ("segmentIndex", ctypes.c_uint32),
                ("triggerIndex", ctypes.c_uint32),
                ("triggerTime", ctypes.c_int64),
                ("timeUnits", ctypes.c_int16),
                ("reserved0", ctypes.c_int16),
                ("timeStampCounter", ctypes.c_uint64)]

  def __init__(self, num_channels, sample_period, raspi=False):
    self.state = 1  # 1 analyse, 2 buy, 3 sell
    self.db = None
    #self.connectDB()
    self.num_channels = num_channels
    self.ps = ps
    self.chandle = ctypes.c_int16()
    self.enabled = 1
    self.disabled = 0
    self.analogue_offset = 0.0
    self.channel_range = self.ps.PS5000A_RANGE['PS5000A_5V']
    self.status = {}
    self.MaxSegments = ctypes.c_uint32()
    #self.buff_size = ctypes.c_int32(buff_size)
    self.noCaptures = ctypes.c_uint32(5)
    self.MaxSamples = ctypes.c_int32()
    #self.triggerInfo = self.PS5000A_TRIGGER_INFO()
    self.timeIntervalns = 0
    self.fromSegmentIndex = ctypes.c_uint32()
    self.toSegmentIndex = ctypes.c_uint32()
    self.timesUpperBits = ctypes.c_uint32()
    self.timesLowerBits = ctypes.c_uint32()
    self.maxADC = 0
    self.sample_period = sample_period
    self.sizeOfOneBuffer = 0
    self.totalSamples = 0
    #self.timebase = 
    self.calculate_time_base()
    self.raspi = False
    if raspi:
      self.init_raspi()
      print("------------raspi set to true")
      self.raspi = True

  def init_raspi(self):
    try:
      #import RPi.GPIO as GPIO
      GPIO_Flag = True
      GPIO.setmode(GPIO.BCM)
      GPIO.setup(21,GPIO.OUT) 
      GPIO.output(21,GPIO.LOW)
      while GPIO.input(21) == 1:
        print("waiting for low")
      print("Raspberry pi Pin 21 enabled and LOW")
    except:
      GPIO_Flag = False
      print("Not raspberryPi GPIO disabled")

  def calculate_time_base(self):
    self.timebase = ctypes.c_uint32(int((self.sample_period * 62.5e6)  + 3)) # from programmers guide formula 5000 series pg 28
  
  def cursorOpen(self):
    try:
      return self.db.cursorOpen()
    except:
      return False

  def time_units(self,sampleTimeInNS):
    # Begin streaming mode:
    self.sampleInterval = ctypes.c_int32(sampleTimeInNS)
    self.sampleUnits = ps.PS5000A_TIME_UNITS['PS5000A_NS']

  def run_streaming(self):
    maxPreTriggerSamples = 0
    autoStopOn = 0
    # No downsampling:
    status = {}
    downsampleRatio = 1
    status["runStreaming"] = self.ps.ps5000aRunStreaming(self.chandle,
                                                    ctypes.byref(self.sampleInterval),
                                                    self.sampleUnits,
                                                    maxPreTriggerSamples,
                                                    0,#self.totalSamples,#totalSamples, Not needed in continous streaming mode - used for autoStop
                                                    autoStopOn,
                                                    downsampleRatio,
                                                    self.ps.PS5000A_RATIO_MODE['PS5000A_RATIO_MODE_NONE'],
                                                    self.sizeOfOneBuffer)
    assert_pico_ok(status["runStreaming"])

  def get_latests_stream_values(self, cFuncPtr):
    status = {}
    status["getStreamingLastestValues"] = ps.ps5000aGetStreamingLatestValues(self.chandle, cFuncPtr, None)
    assert_pico_ok(status["getStreamingLastestValues"])

  def fireTrigger(self, output_pin, delay, pulses, pulse_width_, pulse_width_decay, start_time):
    
    while True:
      if time.time() > (start_time + delay):
        """ GPIO.output(output_pin,GPIO.HIGH)
        while(time.time() < (start_time + delay + 1)):
          #GPIO.output(output_pin,GPIO.HIGH) """
        print(f"Trigger time pulse @@@ ::: {time.time()}")
        break


  def picoSetupHandle(self):
    # Open PicoScope 2000 Series device
    # Resolution set to 12 Bit
    print("Setting up Pico handle...")
    resolution = self.ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_12BIT"]
    # Returns handle to chandle for use in future API functions
    self.status["openunit"] = self.ps.ps5000aOpenUnit(ctypes.byref(self.chandle), None, resolution)

    try:
        assert_pico_ok(self.status["openunit"])
    except: # PicoNotOkError:

        powerStatus = self.status["openunit"]

        if powerStatus == 286:
            self.status["changePowerSource"] = self.ps.ps5000aChangePowerSource(self.chandle, powerStatus)
        elif powerStatus == 282:
            self.status["changePowerSource"] = self.ps.ps5000aChangePowerSource(self.chandle, powerStatus)
        else:
            raise

        assert_pico_ok(self.status["changePowerSource"])
    print("Success")

  def stop_and_close(self):
    status = {}
    status["stop"] = self.ps.ps5000aStop(self.chandle)
    assert_pico_ok(status["stop"])
    print("Successfully executed")
    # Disconnect the scope
    # handle = chandle
    status["close"] = self.ps.ps5000aCloseUnit(self.chandle)
    assert_pico_ok(status["close"])
    if self.raspi:
      GPIO.cleanup()
  
  def picoSetupChannel(self, channel):
    # Set up channel A
    # handle = chandle
    # channel = PS5000A_CHANNEL_A = 0
    # enabled = 1
    # coupling type = PS5000A_DC = 1
    # range = PS5000A_2V = 7
    
    # analogue offset = 0 V
    status = {}
    print(f"Setting up pico channel {channel}")
    #channel_range = self.ps.PS5000A_RANGE['PS5000A_5V']
    status["setChA"] = self.ps.ps5000aSetChannel(self.chandle,
                                            self.ps.PS5000A_CHANNEL[f'PS5000A_CHANNEL_{channel}'],
                                            self.enabled,
                                            self.ps.PS5000A_COUPLING['PS5000A_DC'],
                                            self.channel_range,
                                            self.analogue_offset)
    #print(status)
    assert_pico_ok(status["setChA"])

  def getMaxSegments(self):
    status = {}
    print(f"Getting max mem segments")
    # ctypes.byref(self.MaxSegments)
    status["GetMaxSegments"] = self.ps.ps5000aGetMaxSegments(self.chandle, ctypes.byref(self.MaxSegments))
    print(f"max segments allowed: {self.MaxSegments}")
    assert_pico_ok(status["GetMaxSegments"])

  def setTimeBase(self, buffer_lenght):
    status = {}
    #maxsamples = preTriggerSamples + postTriggerSamples
    # Gets timebase innfomation
    # Handle = chandle
    #self.timebase = n
    # Nosample = maxsamples
    # TimeIntervalNanoseconds = ctypes.byref(timeIntervalns)
    # MaxSamples = ctypes.byref(returnedMaxSamples)
    # Segement index = 0
    print(f"Time base (n): {self.timebase}")
    self.timeIntervalns = ctypes.c_float()
    self.returnedMaxSamples = ctypes.c_int16()

    status["GetTimebase2"] = self.ps.ps5000aGetTimebase2(  self.chandle, 
                                                          self.timebase, 
                                                          buffer_lenght, 
                                                          ctypes.byref(self.timeIntervalns), 
                                                          ctypes.byref(self.returnedMaxSamples), 
                                                          0
                                                        )
    assert_pico_ok(status["GetTimebase2"])




  def getTriggerInfoBulk(self):
    status = {}
    print(f"Getting trigger info bulk")
    # ctypes.byref(self.MaxSegments)
    status["GetTriggerInfoBulk"] = self.ps.ps5000aGetTriggerInfoBulk( self.chandle, 
                                                                      ctypes.byref(self.triggerInfo),   
                                                                      self.fromSegmentIndex, 
                                                                      self.toSegmentIndex)
    #print(f"max segments allowed: {self.MaxSegments}")
    assert_pico_ok(status["GetTriggerInfoBulk"])

  def MemorySegments(self):
    status = {}
    print(f"Setting memory segments")
    # ctypes.byref(self.MaxSegments)
    status["GetMaxSegments"] = self.ps.ps5000aMemorySegments(
                                                              self.chandle, 
                                                              self.MaxSegments, 
                                                              ctypes.byref(self.MaxSamples)
                                                            )
    print(f"max samples allowed: {self.MaxSamples}")
    assert_pico_ok(status["GetMaxSegments"])

  def setNoCaptures(self):
    status = {}
    print(f"Setting max no of captures")
    # ctypes.byref(self.MaxSegments)
    status["SetNoOfCaptures"] = self.ps.ps5000aSetNoOfCaptures(self.chandle, self.noCaptures)
    print(f"SetNoOfCaptures : {self.noCaptures}")
    assert_pico_ok(status["SetNoOfCaptures"])

  def setTrigger(self, delay):
    status = {}
    self.maxADC = ctypes.c_int16()
    status["maximumValue"] = self.ps.ps5000aMaximumValue(self.chandle, ctypes.byref(self.maxADC))
    threshold = int(mV2adc(1000, self.channel_range, self.maxADC))
    print(f"Threhold : {threshold}")
    PS5000A_RISING = 2 # enum is 2 rsiing or falling
    trigger_timeout = 10000 # never
    # Set up single trigger
    # handle = chandle
    # enabled = 1
    # source = PS4000a_CHANNEL_A = 0
    # threshold = 1024 ADC counts
    # direction = PS4000a_RISING = 2
    # delay = 0 s
    # auto Trigger = 1000 ms
    #ps.ps4000aSetSimpleTrigger(chandle, 1, 0, 1024, 2, 0, 100)
    status["SetSimpleTrigger"] = self.ps.ps5000aSetSimpleTrigger(
                                                                            self.chandle, 
                                                                            1, 
                                                                            self.ps.PS5000A_CHANNEL[f'PS5000A_CHANNEL_D'], 
                                                                            threshold, 
                                                                            PS5000A_RISING, 
                                                                            delay, 
                                                                            trigger_timeout)
    #print(f"SetSimpleTrigger set")
    
    assert_pico_ok(status["SetSimpleTrigger"])

  def runStreaming(self, sampleInterval, sampleUnits, totalSamples, sizeOfOneBuffer):
    maxPreTriggerSamples = 1
    autoStopOn = 0
    # No downsampling:
    downsampleRatio = 1
    #status = {}
    self.status["runStreaming"] = self.ps.ps5000aRunStreaming(self.chandle,
                                                    ctypes.byref(sampleInterval),
                                                    sampleUnits,
                                                    maxPreTriggerSamples,
                                                    sizeOfOneBuffer,#totalSamples,
                                                    autoStopOn,
                                                    downsampleRatio,
                                                    self.ps.PS5000A_RATIO_MODE['PS5000A_RATIO_MODE_NONE'],
                                                    sizeOfOneBuffer)
    assert_pico_ok(self.status["runStreaming"])

  def picoSetupChannelBuffer(self, channel, buffer):
    memory_segment = 0
    # Set data buffer location for data collection from channel A
    # handle = chandle
    # source = PS5000A_CHANNEL_A = 0
    # pointer to buffer max = ctypes.byref(bufferAMax)
    # pointer to buffer min = ctypes.byref(bufferAMin)
    # buffer length = maxSamples
    # segment index = 0
    # ratio mode = PS5000A_RATIO_MODE_NONE = 0
    self.status[f"setDataBuffers{channel}"] = self.ps.ps5000aSetDataBuffer(self.chandle,
                                                        self.ps.PS5000A_CHANNEL[f'PS5000A_CHANNEL_{channel}'],
                                                        buffer.ctypes.data_as(ctypes.POINTER(ctypes.c_int16)),
                                                        self.sizeOfOneBuffer,
                                                        memory_segment,
                                                        self.ps.PS5000A_RATIO_MODE['PS5000A_RATIO_MODE_NONE'])
    assert_pico_ok(self.status[f"setDataBuffers{channel}"])

  def getValuesTriggerTimeBulk(self):
    memory_segment = 0
    # Set data buffer location for data collection from channel A
    # handle = chandle
    # source = PS5000A_CHANNEL_A = 0
    # pointer to buffer max = ctypes.byref(bufferAMax)
    # pointer to buffer min = ctypes.byref(bufferAMin)
    # buffer length = maxSamples
    # segment index = 0
    # ratio mode = PS5000A_RATIO_MODE_NONE = 0
    self.status[f"getValuesTriggerTimeBulk"] = self.ps.ps5000aGetValuesTriggerTimeOffsetBulk(  self.chandle,
                                                                                              ctypes.byref(self.timesUpperBits),
                                                                                              ctypes.byref(self.timesUpperBits),
                                                                                              self.ps.PS5000A_TIME_UNITS['PS5000A_NS'],
                                                                                              self.fromSegmentIndex,
                                                                                              self.toSegmentIndex)
    assert_pico_ok(self.status[f"getValuesTriggerTimeBulk"])
      
  def connectDB(self):
    try:
      time.sleep(1)
      db_name = "pi_measure" 
      ip = "127.0.0.1"
      user = "root"
      password = "password"
      self.db = database.database(db_name, user, password, ip)
    except Exception as err:
      print(f"Error connecting to database. Error: {err}")
      raise err


  def timeWindow(self, stamp, delay):
    delta = 0
    while delta < delay:
      time.sleep(1)
      delta = (datetime.datetime.now() - stamp).seconds

  def getLogs(self, table_name, time="1 HOUR", reverse=True, raw=False):
    data, code = self.db.getAllRows(table_name, time=time, reverse=reverse, raw=raw)
    return data

  def format_call(self, name, results):
    print('-'*80)
    print('%50s' % (name,))
    print('-'*80)
    pp.pprint(results)
    print('-'*80)

