from numpy import *
from scipy.integrate import cumtrapz
import time
import math

"""v11 - has additional input to snakecall that means the GUI v2.5 can call
custom angles and velocityies for specific joint control"""

def x_cart(elevation,azimuth,r):                   #return the x coordinates
  return cos(elevation)*cos(azimuth)*(r)

def y_cart(elevation,azimuth,r):                   #return the y coordinates
   return cos(elevation)*sin(azimuth)*(r)

def z_cart(elevation,r):                           #return the z coordinates
   return sin(elevation)*(r)

def sigmf(x,a,c):
    return 1./(1 + exp(-a*(x-c)));

def decodeTextFile(Filename):
   """This function takes a text file with two columns of 10 numbers and decodes it
      into a 2d array which can be interpreted as angles"""
   angleList = list()
   f = open(Filename,"r")
   fText = f.read()
   angles = fText.split("\n")
   print angles
   for i in range(0,10):
      angleList.append(angles[i].split(" "))
   print angleList
   return angleList
   f.close()

def writeCustomAngles(ang_and_vel_List,instruction):
   AnglesFile = open("snakeAngles.txt",instruction)
   VelsFile = open("snakeVelocities.txt",instruction)
   print "called"
   for section in range(0,len(ang_and_vel_List)):
     for val in range(0,4):
       if ang_and_vel_List[section][val]=="":
         ang_and_vel_List[section][val]="0.0"
     AnglesFile.write("%s %s\n" %(ang_and_vel_List[section][0],ang_and_vel_List[section][1]))
     VelsFile.write("%s %s\n" %(ang_and_vel_List[section][2],ang_and_vel_List[section][3]))
   AnglesFile.close()
   VelsFile.close()
        
def snake(total_time,Ahor,Aver,rhohor,rhover,deltahor,deltaver,betahor,betaver,delta0,reset=False,instruction='w'):
  """return thetahor and thetaver"""                               
  
  N = 11                              #number of modules
  S = 0.25                             #length of each modules (20cm)
  dt = 1.0                            #time step length
  
  steps=floor(total_time/dt)          #number of time steps - rounded to
  total_time = steps*dt
  thetaver = zeros((steps,N))         #rotation angle about vertical axis
  thetahor = zeros((steps,N))         #rotation angle about horizontal axis

  nn = 100                            # number of points that describes the S-curves
  ang_vhor = zeros((steps,N-1))
  ang_vver = zeros((steps,N-1))
  
  ang_vhorn = zeros((steps,N-1,nn))
  ang_vvern = zeros((steps,N-1,nn))

  ang_hor = zeros((steps,N-1,nn))
  ang_ver = zeros((steps,N-1,nn))
  
  Tx = zeros((steps,N-1,nn))
  Ty = zeros((steps,N-1,nn))
  tn = linspace(0,dt,nn)
  sigmoid = array((sigmf(tn,15,dt/2))) #makes the s-curve
  x = zeros((steps,N+1))              #X coordinate of joints
  y = zeros((steps,N+1))              #Y coordinate of joints
  z = zeros((steps,N+1))              #Z coordinate of joints
  CoGx = zeros((steps,N))             #x coordinate of CoG
  CoGy = zeros((steps,N))             #y coordinate of CoG
  f=open('snakeAngles.txt',instruction)            #Open/create a text file for exporting thetahor thetaver
  g=open('snakeVelocities.txt',instruction)        #Open/create a text file for exporting angvelohor angvelover
  h=open('snakeTorques.txt',instruction)
  initialAngles =decodeTextFile("SnakeAngleFeedback.txt")
  n=0                                 #counter of time steps
  
  m = 5                               # mass of each segment (kg)
  I = (m*S**2)/12                     # moment of inertia (for a thin rod)
  
  reset=False
  for t in linspace(dt,total_time,steps):
    "Changed to go form dt - as in Leonard's code"
    phi=0                             #gravity vector can change with time
    for i in range(N):
      if reset == False:
        horAng = betahor+Ahor*sin(rhohor*t+i*deltahor+delta0)
        verAng = betaver+Aver*sin(rhover*t+i*deltaver)
        thetahor[n,i]=cos(phi)*horAng - sin(phi)*verAng       #calculate thetahor
        thetaver[n,i]=sin(phi)*horAng + cos(phi)*verAng       #calculate thetaver
      else:
        thetahor[n,i]=0
        thetaver[n,i]=0
      
      x[n,i+1] = x_cart(thetaver[n,i], thetahor[n,i],S) + x[n,i]      #calculate X coordinate of joints
      y[n,i+1] = y_cart(thetaver[n,i], thetahor[n,i],S) + y[n,i]      #calculate Y coordinate of joints
      z[n,i+1] = z_cart(thetaver[n,i],S) + z[n,i]                     #calculate Z coordinate of joints
      CoGx[n,i] = (x[n,i+1]+x[n,i])/2                                 #calculate X coordinate of COG of joints
      CoGy[n,i] = (y[n,i+1]+y[n,i])/2                                 #calculate Y coordinate of COG of joints

    for i in range(N-1):
      f.write("% f " %thetaver[n,i])
      f.write("% f\n" %thetahor[n,i])
      if reset==False:
        if n == 0: #first iteration -> compare angles to current state of snake
          #+ feedback rather than - because of polarity difference between sim and ME
          # 0th column in feedback data is VERTICAL data, 1st is horizontal.
          if (abs(float(initialAngles[i][0]))+abs(float(initialAngles[i][0])))>0:
            ang_vhor[n,i]=0.4
            ang_vver[n,i]=0.4
          else:
            ang_vhor[n,i]=(thetahor[n,i]+float(initialAngles[i][1]))/dt
            print "theta_hor = %f; fb = %f; vel = %f" %(thetahor[n,i],-float(initialAngles[i][1]),ang_vhor[n,i])
            ang_vver[n,i]=(thetaver[n,i]+float(initialAngles[i][0]))/dt
            print "theta_ver = %f; fb = %f; vel = %f" %(thetaver[n,i],-float(initialAngles[i][0]),ang_vver[n,i])
        else:
          ang_vhor[n,i]=(thetahor[n,i]-thetahor[n-1,i])/dt
          ang_vver[n,i]=(thetaver[n,i]-thetaver[n-1,i])/dt
        g.write("% f " %ang_vver[n,i])
        g.write("% f\n" %ang_vhor[n,i])
      else:
        g.write("% f " % (1))        # reset with linear ang vel of 1.0   
        g.write("% f\n" % (1))

      max_horvel = 0.5*(thetahor[n,i] + thetahor[n-1,i])
      max_vervel = 0.5*(thetaver[n,i] + thetaver[n-1,i])
      
      sub_points = 20
      
      "+ve grad vel curve"
      
      
      

      """n = number of time steps (int), i = joint number"""   
      ang_vhorn[n,i,:]=ang_vhor[n,i]*sigmoid
      ang_vvern[n,i,:]=ang_vver[n,i]*sigmoid
      """
      ang_hor[n,i,:]=(append([0],cumtrapz(ang_vhorn[n,i,:])))+thetahor[n,i] #just needed for inbetween angle data
      ang_ver[n,i,:]=(append([0],cumtrapz(ang_vvern[n,i,:])))+thetaver[n,i] #feel free to ignore both, doesn't do anything practical
      """
      ang_hor[n,i,:]=thetahor[n,i]*sigmoid + thetahor[n-1,i] #offset terms from previous angle calculated
      ang_ver[n,i,:]=thetaver[n,i]*sigmoid + thetaver[n-1,i] #future work - should be angle feedback instead of offset
      
      alphahor=diff(ang_vhorn[n,i,:]/(dt/nn))
      alphaver=diff(ang_vvern[n,i,:]/(dt/nn))
      alphahor=append(alphahor,-alphahor[0])
      alphaver=append(alphaver,-alphaver[0])
      Tx[n,i,:]=sum((I + m*(CoGx[n,:]-x[n,i])**2))*alphahor         # torque around joint, taking anti-clockwise as +ve
      Ty[n,i,:]=sum((I + m*(CoGy[n,:]-y[n,i])**2))*alphaver
      h.write("% f " % Tx[n,i,-1])        #put : on the nn index to get all the data   
      h.write("% f\n" % Ty[n,i,-1])
    n=n+1


  f.close()
  g.close()
  h.close()
  
############################################################################################################################################
def halfsnake(total_time,AhorF,AhorR,AverF,AverR,rhohorF,rhohorR,rhoverF,rhoverR,deltahorF,deltahorR,deltaverF,deltaverR,betahorF,betahorR,betaverF,betaverR,delta0,reset=False,instruction='w'):
  """return thetahor and thetaver (horizontal and vertical rotational angles)where
            Ahor=horizontal wave amplitude,Aver=vertical wave amplitude,
            rhohor=horizontal wave angular frequency,rhover= vertical wave angular frequency,
            deltahor=horzontal wave phase difference bettween modules,
            deltaver=vertical wave phase difference bettween modules,
            betahor=horizontal wave phase offset,
            betaver=vertical wave phase offset,
            delta0=phase offset between ver & hor waves
            reset = whether you want to return 0's for all angles
            instruction = to re-write the text file or append to it (default='a')"""                                
  
  N = 11                              #number of modules
  S = 0.25                             #length of each modules (20cm)
  dt = 1.0                            #time step length
  steps=total_time/dt                 #number of time steps
  verAng = zeros((steps,N))         #rotation angle about vertical axis
  horAng = zeros((steps,N))         #rotation angle about horizontal axis
  thetaver = zeros((steps,N))         #rotation angle about vertical axis
  thetahor = zeros((steps,N))         #rotation angle about horizontal axis

  nn = 100                            # number of points that describes the S-curves
  ang_vhor = zeros((steps,N-1))
  ang_vver = zeros((steps,N-1))
  
  ang_vhorn = zeros((steps,N-1,nn))
  ang_vvern = zeros((steps,N-1,nn))

  ang_hor = zeros((steps,N-1,nn))
  ang_ver = zeros((steps,N-1,nn))
  
  Tx = zeros((steps,N-1,nn))
  Ty = zeros((steps,N-1,nn))
  tn = linspace(0,dt,nn)
  sigmoid = array((sigmf(tn,15,dt/2))) #makes the s-curve
  x = zeros((steps,N+1))              #X coordinate of joints
  y = zeros((steps,N+1))              #Y coordinate of joints
  z = zeros((steps,N+1))              #Z coordinate of joints
  CoGx = zeros((steps,N))             #x coordinate of CoG
  CoGy = zeros((steps,N))             #y coordinate of CoG
  f=open('snakeAngles.txt',instruction)            #Open/create a text file for exporting thetahor thetaver
  g=open('snakeVelocities.txt',instruction)        #Open/create a text file for exporting angvelohor angvelover
  h=open('snakeTorques.txt',instruction)
  initialAngles =decodeTextFile("SnakeAngleFeedback.txt")
  n=0                                 #counter of time steps
  
  m = 5                               # mass of each segment (kg)
  I = (m*S**2)/12                     # moment of inertia (for a thin rod)
  
  reset=False
  for t in linspace(dt,total_time,steps):
    "Changed to go form dt - as in Leonard's code"
    phi=0                             #gravity vector can change with time
    for i in range(N-6):
      if reset == False:
        horAng[n,i] = betahorF+AhorF*sin(rhohorF*t+i*deltahorF+delta0)
        verAng[n,i] = betaverF+AverF*sin(rhoverF*t+i*deltaverF)
        thetahor[n,i]=cos(phi)*horAng[n,i] - sin(phi)*verAng[n,i]      #calculate thetahor
        thetaver[n,i]=sin(phi)*horAng[n,i]+ cos(phi)*verAng[n,i]      #calculate thetaver
      else:
        thetahor[n,i]=0
        thetaver[n,i]=0

    for i in range(N-6,N):
      if reset == False:
        horAng[n,i] = betahorR+AhorR*sin(rhohorR*t+i*deltahorR+delta0)
        verAng[n,i]= betaverR+AverR*sin(rhoverR*t+i*deltaverR)
        thetahor[n,i]=cos(phi)*horAng[n,i] - sin(phi)*verAng[n,i]       #calculate thetahor
        thetaver[n,i]=sin(phi)*horAng[n,i] + cos(phi)*verAng[n,i]       #calculate thetaver
      else:
        thetahor[n,i]=0
        thetaver[n,i]=0
        
      x[n,i+1] = x_cart(thetaver[n,i], thetahor[n,i],S) + x[n,i]      #calculate X coordinate of joints
      y[n,i+1] = y_cart(thetaver[n,i], thetahor[n,i],S) + y[n,i]      #calculate Y coordinate of joints
      z[n,i+1] = z_cart(thetaver[n,i],S) + z[n,i]                     #calculate Z coordinate of joints
      CoGx[n,i] = (x[n,i+1]+x[n,i])/2                                 #calculate X coordinate of COG of joints
      CoGy[n,i] = (y[n,i+1]+y[n,i])/2                                 #calculate Y coordinate of COG of joints

    for i in range(N-1):
      f.write("% f " %thetaver[n,i])
      f.write("% f\n" %thetahor[n,i])
      if reset==False:
        if n == 0: #first iteration - compare angles to current state of snake
          #+ve float because of inverse polarity with simulation in ME.
          #note index of file being read also different as [0] are vertical angles.
          ang_vhor[n,i]=(thetahor[n,i]+float(initialAngles[i][1]))/dt
          ang_vver[n,i]=(thetaver[n,i]+float(initialAngles[i][0]))/dt
        else:
          ang_vhor[n,i]=(thetahor[n,i]-thetahor[n-1,i])/dt
          ang_vver[n,i]=(thetaver[n,i]-thetaver[n-1,i])/dt
        g.write("% f " %ang_vver[n,i])
        g.write("% f\n" %ang_vhor[n,i])
      else:
        g.write("% f " % (1))        # reset with linear ang vel of 1.0   
        g.write("% f\n" % (1))
      
      max_horvel = 0.5*(thetahor[n,i] + thetahor[n-1,i])
      max_vervel = 0.5*(thetaver[n,i] + thetaver[n-1,i])
      
      sub_points = 20
      
      "+ve grad vel curve"
      
      
      

      """n = number of time steps (int), i = joint number"""   
      ang_vhorn[n,i,:]=ang_vhor[n,i]*sigmoid
      ang_vvern[n,i,:]=ang_vver[n,i]*sigmoid
      """
      ang_hor[n,i,:]=(append([0],cumtrapz(ang_vhorn[n,i,:])))+thetahor[n,i] #just needed for inbetween angle data
      ang_ver[n,i,:]=(append([0],cumtrapz(ang_vvern[n,i,:])))+thetaver[n,i] #feel free to ignore both, doesn't do anything practical
      """
      ang_hor[n,i,:]=thetahor[n,i]*sigmoid + thetahor[n-1,i] #offset terms from previous angle calculated
      ang_ver[n,i,:]=thetaver[n,i]*sigmoid + thetaver[n-1,i] #future work - should be angle feedback instead of offset
      
      alphahor=diff(ang_vhorn[n,i,:]/(dt/nn))
      alphaver=diff(ang_vvern[n,i,:]/(dt/nn))
      alphahor=append(alphahor,-alphahor[0])
      alphaver=append(alphaver,-alphaver[0])
      Tx[n,i,:]=sum((I + m*(CoGx[n,:]-x[n,i])**2))*alphahor         # torque around joint, taking anti-clockwise as +ve
      Ty[n,i,:]=sum((I + m*(CoGy[n,:]-y[n,i])**2))*alphaver
      h.write("% f " % Tx[n,i,-1])        #put : on the nn index to get all the data   
      h.write("% f\n" % Ty[n,i,-1])
    n=n+1


  f.close()
  g.close()
  h.close()
  
################################################################################################################################  

def snakecall(gaits,time=20,orientation='0',speed='5',reset=False,instruction="w",custAng=[]):
    """ 'types' = 'sidewinding'/'lateral'/'rotating'/'rolling'/'linear'.
    'forward' means directly sideway for sidewinding where delta0 = pi/4"""
    if len(custAng)==0:
        if reset == True:
           return snake(20,pi/6,pi/6,3*pi/4,3*pi/4,7*pi/18,7*pi/18,0,0,pi/4,True,instruction)
                    
        if gaits == 'sidewinding':
                if orientation == '-1':
                   if speed == '3':
                      return snake(time,pi/6,pi/6,3*pi/4,3*pi/4,29*pi/72,29*pi/72,0,0,71*pi/180,False,instruction)
                   elif speed == '2':
                      return snake(time,pi/6,pi/6,3*pi/4,3*pi/4,4*pi/9,4*pi/9,0,0,11*pi/30,False,instruction)
                   elif speed == '1':
                      return snake(time,pi/6,pi/6,3*pi/4,3*pi/4,17*pi/36,17*pi/36,0,0,61*pi/180,False,instruction)
                   else:
                     print "%s is an invalid speed for the %s gait and %s orientation, please choose a value between 1 and 3" %(speed,gaits,orientation)
                elif orientation == '1':
                   if speed == '3':
                      return snake(time,pi/6,pi/6,3*pi/4,3*pi/4,29*pi/72,29*pi/72,0,0,-71*pi/180,False,instruction)
                   elif speed == '2':
                      return snake(time,pi/6,pi/6,3*pi/4,3*pi/4,4*pi/9,4*pi/9,0,0,-13*pi/36,False,instruction)
                   elif speed == '1':
                      return snake(time,pi/6,pi/6,3*pi/4,3*pi/4,17*pi/36,17*pi/36,0,0,-61*pi/180,False,instruction)
                   else:
                     print "%s is an invalid speed for the %s gait and %s orientation, please choose a value between 1 and 3" %(speed,gaits,orientation) 
                else:
                   print "%s is an invalid orientation for the %s gait, please either '1' or '-1'" %(orientation,gaits)
                    
 

        elif gaits == 'rotating':
                if orientation == '-1':
                   if speed == '0':
                      return halfsnake(time,pi/6,-pi/6,pi/6,-pi/6,3*pi/4,-3*pi/4,3*pi/4,-3*pi/4,29*pi/72,-29*pi/72,29*pi/72,-29*pi/72,0,0,0,0,pi/4,False,instruction)
                   elif speed == '1':
                      return halfsnake(time,pi/6,-pi/6,pi/6,-pi/6,3*pi/4,-3*pi/4,3*pi/4,-3*pi/4,5*pi/16,-5*pi/16,5*pi/16,-5*pi/16,0,0,0,0,pi/4,False,instruction)
                   else:
                     print "%s is an invalid speed for the %s gait and %s orientation, please choose either '0' or '1'" %(speed,gaits,orientation) 
                elif orientation == '1':
                     if speed == '0':
                        return halfsnake(time,-pi/6,pi/6,-pi/6,pi/6,-3*pi/4,3*pi/4,-3*pi/4,3*pi/4,-29*pi/72,29*pi/72,-29*pi/72,29*pi/72,0,0,0,0,3*pi/4,False,instruction)
                     elif speed == '1':
                        return halfsnake(time,-pi/6,pi/6,-pi/6,pi/6,-3*pi/4,3*pi/4,-3*pi/4,3*pi/4,-5*pi/16,5*pi/16,-5*pi/16,5*pi/16,0,0,0,0,3*pi/4,False,instruction)
                     else:
                        print "%s is an invalid speed for the %s gait and %s orientation, please choose either '0' or '1'" %(speed,gaits,orientation)
                else:
                   print "%s is an invalid orientation for the %s gait, please choose either '-1' or '1'" %(orientation,gaits)
              
        elif gaits == 'rolling':
                #orientation = +1 goes right 
                if orientation == '1':
                   if speed == '3':
                      return snake(time,1*pi/18,1*pi/18,5*pi/6,5*pi/6,0,0,0,0,4*pi/9,False,instruction)
                   elif speed == '2':
                      return snake(time,1*pi/18,1*pi/18,13*pi/18,13*pi/18,0,0,0,0,4*pi/9,False,instruction)
                   elif speed == '1':
                      return snake(time,1*pi/18,1*pi/18,10*pi/18,10*pi/18,0,0,0,0,4*pi/9,False,instruction)                    
                   else:
                      print "%s is an invalid speed for the %s gait and %s orientation, please choose a value between 1 and 5" %(speed,gaits,orientation)
                #orientation = -1 goes left
                elif orientation == '-1':
                   if speed == '3':
                      return snake(time,1*pi/18,1*pi/18,5*pi/6,5*pi/6,0,0,0,0,-4*pi/9,False,instruction)
                   elif speed == '2':
                      return snake(time,1*pi/18,1*pi/18,13*pi/18,13*pi/18,0,0,0,0,-4*pi/9,False,instruction)
                   elif speed == '1':
                      return snake(time,1*pi/18,1*pi/18,10*pi/18,10*pi/18,0,0,0,0,-4*pi/9,False,instruction)                
                   else:
                      print "%s is an invalid speed for the %s gait and %s orientation, please choose a value between 1 and 5" %(speed,gaits,orientation)
                else:
                   print "%s is an invalid orientation for the %s gait, please choose either '1' or '-1'" %(orientation,gaits)
              
        elif gaits == 'linear':
                #forwards
                if orientation == '1':
                   if speed == '1':
                      return halfsnake(time,0,0,pi/4,pi/4,0,0,5*pi/6,5*pi/6,0,0,2*pi/3,2*pi/3,-23*pi/360,23*pi/360,0,0,0,False,instruction)   
                   elif speed == '2':
                      return halfsnake(time,0,0,13*pi/48,13*pi/48,0,0,5*pi/6,5*pi/6,0,0,2*pi/3,2*pi/3,-1*pi/15,1*pi/15,0,0,0,False,instruction)
                   elif speed == '3':
                      return halfsnake(time,0,0,7*pi/24,7*pi/24,0,0,5*pi/6,5*pi/6,0,0,2*pi/3,2*pi/3,-1*pi/15,1*pi/15,0,0,0,False,instruction)                    
                   else:
                      print "%s is an invalid speed for the %s gait and %s orientation, please choose a value between 1 and 3" %(speed,gaits,orientation) 
                #backwards
                if orientation == '-1':
                   if speed == '1':                     
                      return halfsnake(time,0,0,pi/4,pi/4,0,0,-5*pi/6,-5*pi/6,0,0,2*pi/3,2*pi/3,23*pi/360,-23*pi/360,0,0,0,False,instruction)
                   elif speed == '2':
                      return halfsnake(time,0,0,13*pi/48,13*pi/48,0,0,-5*pi/6,-5*pi/6,0,0,2*pi/3,2*pi/3,-23*pi/360,23*pi/360,0,0,0,False,instruction)                     
                   elif speed == '3':
                      return halfsnake(time,0,0,7*pi/24,7*pi/24,0,0,-5*pi/6,-5*pi/6,0,0,2*pi/3,2*pi/3,-23*pi/360,23*pi/360,0,0,0,False,instruction)  
                   else:
                      print "%s is an invalid speed for the %s gait and %s orientation, please choose a value between 1 and 3" %(speed,gaits,orientation)
                else:
                   print "%s is an invalid orientation for the %s gait, please choose one of the following:\n'00'\n'01'\n'-01'\n'10'\n'-11'\n'11'" %(orientation,gaits)
                
        elif gaits == 'slithering':
                if orientation== '-1':
                   if speed == '1':                  
                      return snake(time,pi/4,pi/4,5*pi/12,5*pi/6,-5*pi/18,-3*pi/4,0,0,4*pi/9,False,instruction)
                   elif speed == '2':                  
                      return snake(time,5*pi/18,5*pi/18,5*pi/12,5*pi/6,-5*pi/18,-3*pi/4,0,0,19*pi/45,False,instruction)                    
                   elif speed == '3':                  
                      return snake(time,pi/3,pi/3,5*pi/12,5*pi/6,-5*pi/18,-3*pi/4,0,0,13*pi/30,False,instruction)
                   else:
                      print "%s is an invalid speed for the %s gait and %s orientation, please choose a value between 1 and 3" %(speed,gaits,orientation)
                    
                elif orientation== '1':
                   if speed == '1':                  
                      return snake(time,pi/4,pi/4,5*pi/12,5*pi/6,5*pi/18,3*pi/4,0,0,pi/4,False,instruction)
                   elif speed == '2':                  
                      return snake(time,5*pi/18,5*pi/18,5*pi/12,5*pi/6,5*pi/18,3*pi/4,0,0,19*pi/45,False,instruction)                    
                   elif speed == '3':                  
                      return snake(time,pi/3,pi/3,5*pi/12,5*pi/6,5*pi/18,3*pi/4,0,0,13*pi/30,False,instruction) 

                   else:
                      print "%s is an invalid speed for the %s gait and %s orientation, please choose a value between 1 and 3" %(speed,gaits,orientation)
                else:
                   print "%s is an invalid orientation for the %s gait, please choose either '1' or '-1'" %(orientation,gaits)
              
        else:
          print "%s is not a valid gait, please re-run the code with one of the following as the gait:\nsidewinding\nlinear\nrolling\nslithering\nrotating\nlateral" %gaits
    else:
        writeCustomAngles(custAng,instruction)
      
