import os, time
class Setup_PI(object):

  def __init__(self, time_stamp):
    #self.time_stamp = time_stamp
    self.sync_clock_pi(time_stamp)

  def sync_clock_pi(self, time_stamp):
      # NTP must be inactive on PI!!! - sudo timedatectl set-ntp false
      #ref_time = float(time_stamp)
      """ st = time.ctime(ref_time)
      list_st = st.split(" ")
      year = list_st.pop()
      list_st.append("CET")
      new_date = ""
      for item in list_st:
        if item == "":
          continue
        new_date += item + " "
      new_date += year
      command = f"sudo date -s '{new_date}'" """
      #command = f"sudo date +%s -s @{time_stamp}"
      #os.system(command)
      #print(f"Date altered successfully -  \t Clock: {os.system('date')}")
      #print(command)
      pass
      #start_object = ref_time
      #rpi = True

  def setup_GPIO(self, output_pin):
    try:
      import RPi.GPIO as GPIO
      GPIO_Flag = True
      GPIO.setmode(GPIO.BCM)
      GPIO.setup(output_pin, GPIO.OUT) 
      GPIO.output(output_pin, GPIO.LOW)
      while GPIO.input(output_pin) == 1:
        print("waiting for low")
      print("Raspberry pi Pin 21 enabled and LOW")
    except:
      GPIO_Flag = False
      print("Not raspberryPi GPIO disabled")