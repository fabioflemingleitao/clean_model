# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 23:40:48 2025

@author: fabio
"""
'''
# Hello everyone,
Sorry for resurrecting this topic.  (https://www.facebook.com/share/p/1DuNexzfdK/) https://www.facebook.com/groups/opensees
However, I haven’t been able to find a solution. I’ve been running analyses 
on frame structures subjected to seismic accelerations. Some analyses converge, 
and others do not (which is normal). However, in some specific situations 
(usually with high acceleration signals due to a high scale factor), a fatal error 
occurs and the analysis is automatically interrupted (crashing the analysis).
I’ve been trying every day to solve this error, but without success. 


By debugging this code, the "fatal error" occurs in an "ops.analyze" command 
around time = 42.421s of the earthquake. Perhaps by debugging inside OpenSees (in C++) 
it might be possible to identify what’s happening. But I don’t know how to do that.

# Today, I decided to create a simple model where I could replicate this "fatal error" 
so you can try to help me.

Note: I noticed that in this analysis, when I remove the reinforcement from my 
beam section, the fatal error does not occur.


Thank you for the help.
Fábio Leitão
'''




import openseespy.opensees as ops
import opsvis as opsv
import numpy as np
import matplotlib.pyplot as plt
import ReadRecord
import openseespy.preprocessing.DiscretizeMember as opsdm



# ops.logFile("novo_log2.txt", 1)
print("Starting RCFrame Clean Model")
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)



# =============================================================================
# 1 - Create nodes
# =============================================================================

width = 4.5
height = 3


ops.node(1, 0, 0)
ops.node(2, width, 0.0)
ops.node(3, 2*width, 0.0)
ops.node(4, 3*width, 0.0)

ops.node(5, 0, height)
ops.node(6, width, height)
ops.node(7, 2*width, height)
ops.node(8, 3*width, height)

ops.node(9, 0, 2*height)
ops.node(10, width, 2*height)
ops.node(11, 2*width, 2*height)
ops.node(12, 3*width, 2*height)

ops.node(13, 0, 3*height)
ops.node(14, width, 3*height)
ops.node(15, 2*width, 3*height)
ops.node(16, 3*width, 3*height)


# =============================================================================
#  2 - Fix nodes base
# =============================================================================

ops.fix(1, 1, 1, 1)
ops.fix(2, 1, 1, 1)
ops.fix(3, 1, 1, 1)
ops.fix(4, 1, 1, 1)


# =============================================================================
#  3 - Define materials for beam and columns
# =============================================================================

# Concrete Non-Confined
ops.uniaxialMaterial("Concrete07", 1, -33009999.999999996, -0.0020795618440196985, 30430278285.47081, 3562168.4407113595, 0.000234120004246701, 2, 2.3, 4.448076923076922)

# Core concrete (confined) - Beam
ops.uniaxialMaterial("Concrete07", 5 ,  -33946950.00079997 ,  -0.002374691421635245 ,  30430278285.47081 ,  3562168.4407113595 ,  0.000234120004246701 ,  2 ,  30 ,  1.885983147073742)

# Core concrete (confined) - Column
ops.uniaxialMaterial("Concrete07", 4 ,  -35041489.76597064 ,  -0.002719460148144062 ,  30430278285.47081 ,  3562168.4407113595 ,  0.000234120004246701 ,  2 ,  30 ,  1.7344309001863512)

# Steel
ops.uniaxialMaterial('Steel02', 2, 576000000, 201000000000, 0.012, 20, 0.9, 0.08, 0.039, 1, 0.029, 1, 0.0)


# =============================================================================
# 4.1 - Define Beam FiberSection
# =============================================================================

b_beam = 0.15
h_beam = 0.40
Ec = 30430278285.47081
cover = 0.04
As_Lower = 3.1415/4 * (0.0125)**2
As_Upper = 3.1415/4 * (0.010)**2

Area = b_beam*h_beam
Iz = b_beam*h_beam**3/12
Iy = h_beam*b_beam**3/12
J = Iz + Iy  
Gc = Ec/(2*(1+0.2))
GJ = Gc * J
   
h = h_beam
b = b_beam
c = cover
NumDivCore_z = int(np.floor((b-2*c)/c))
NumDivCore_y = int(np.floor((h-2*c)/c))
NumDivRight_z = 1
NumDivRight_y = int(np.floor(h/c))
NumDivLeft_z = 1
NumDivLeft_y = int(np.floor(h/c))
NumDivLower_z = int(np.floor(b/c))
NumDivLower_y = 1   
NumDivUpper_z = int(np.floor(b/c))
NumDivUpper_y = 1

fib_sec_1 = [['section', 'Fiber', 1, '-GJ', GJ],
            ['patch', 'quad', 5, NumDivCore_y, NumDivCore_z, -h/2+c, -b/2+c, h/2-c, -b/2+c, h/2-c, b/2-c, -h/2+c, b/2-c],  # noqa: E501
          ['patch', 'quad', 1, NumDivRight_y, NumDivRight_z, -h/2, -b/2, h/2, -b/2, h/2, -b/2+c, -h/2, -b/2+c],
          ['patch', 'quad', 1, NumDivLeft_y, NumDivLeft_z, -h/2, b/2-c, h/2, b/2-c, h/2, b/2, -h/2, b/2],  
          ['patch', 'quad', 1, NumDivLower_y, NumDivLower_z, -h/2, -b/2+c, -h/2+c, -b/2+c, -h/2+c, b/2-c, -h/2,b/2-c],
          ['patch', 'quad', 1, NumDivUpper_y, NumDivUpper_z, h/2-c, -b/2+c, h/2, -b/2+c, h/2, b/2-c, h/2-c, b/2-c],
          ]


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# If you comment out this part (and remove the reinforcement from the beam), the 
# fatal error does not occur for this analysis.

fib_sec_1.append(['layer', 'straight', 2, 2, As_Upper, h/2-c, b/2-c, h/2-c, c-b/2])
fib_sec_1.append(['layer', 'straight', 2, 3, As_Lower, c-h/2, b/2-c, c-h/2, c-b/2])

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

opsv.fib_sec_list_to_cmds(fib_sec_1)
# opsv.plot_fiber_section(fib_sec_1)
# plt.axis('equal')


# =============================================================================
# 4.2 - Define Column Fiber Section
# =============================================================================

b_col = 0.25
h_col = 0.25
Ec = 30430278285.47081
cover = 0.04
As = 3.1415/4 * (0.01)**2

Area = b_col*h_col
Iz = b_col*h_col**3/12
Iy = h_col*b_col**3/12
J = Iz + Iy  
Gc = Ec/(2*(1+0.2))
GJ = Gc * J
   
h = h_col
b = b_col
c = cover
NumDivCore_z = int(np.floor((b-2*c)/c))
NumDivCore_y = int(np.floor((h-2*c)/c))
NumDivRight_z = 1
NumDivRight_y = int(np.floor(h/c))
NumDivLeft_z = 1
NumDivLeft_y = int(np.floor(h/c))
NumDivLower_z = int(np.floor(b/c))
NumDivLower_y = 1   
NumDivUpper_z = int(np.floor(b/c))
NumDivUpper_y = 1

fib_sec_2 = [['section', 'Fiber', 2, '-GJ', GJ],
            ['patch', 'quad', 4, NumDivCore_y, NumDivCore_z, -h/2+c, -b/2+c, h/2-c, -b/2+c, h/2-c, b/2-c, -h/2+c, b/2-c],  # noqa: E501
          ['patch', 'quad', 1, NumDivRight_y, NumDivRight_z, -h/2, -b/2, h/2, -b/2, h/2, -b/2+c, -h/2, -b/2+c],
          ['patch', 'quad', 1, NumDivLeft_y, NumDivLeft_z, -h/2, b/2-c, h/2, b/2-c, h/2, b/2, -h/2, b/2],  
          ['patch', 'quad', 1, NumDivLower_y, NumDivLower_z, -h/2, -b/2+c, -h/2+c, -b/2+c, -h/2+c, b/2-c, -h/2,b/2-c],
          ['patch', 'quad', 1, NumDivUpper_y, NumDivUpper_z, h/2-c, -b/2+c, h/2, -b/2+c, h/2, b/2-c, h/2-c, b/2-c],
          ]

fib_sec_2.append(['layer', 'straight', 2, 2, As, h/2-c, b/2-c, h/2-c, c-b/2])
fib_sec_2.append(['layer', 'straight', 2, 2, As, c-h/2, b/2-c, c-h/2, c-b/2])

opsv.fib_sec_list_to_cmds(fib_sec_2)


# =============================================================================
# 5.1 - Create Beams
# =============================================================================

num_discretize = 5
ops.beamIntegration('Lobatto', 1, 1, 3)
ops.geomTransf('Linear', 1)

opsdm.DiscretizeMember(5, 6, num_discretize, 'dispBeamColumn', 1, 1, 50, 50)
opsdm.DiscretizeMember(6, 7, num_discretize, 'dispBeamColumn', 1, 1, 60, 60)
opsdm.DiscretizeMember(7, 8, num_discretize, 'dispBeamColumn', 1, 1, 70, 70)
opsdm.DiscretizeMember(9, 10, num_discretize, 'dispBeamColumn', 1, 1, 90, 90)
opsdm.DiscretizeMember(10, 11, num_discretize, 'dispBeamColumn', 1, 1, 100, 100)
opsdm.DiscretizeMember(11, 12, num_discretize, 'dispBeamColumn', 1, 1, 110, 110)
opsdm.DiscretizeMember(13, 14, num_discretize, 'dispBeamColumn', 1, 1, 130, 130)
opsdm.DiscretizeMember(14, 15, num_discretize, 'dispBeamColumn', 1, 1, 140, 140)
opsdm.DiscretizeMember(15, 16, num_discretize, 'dispBeamColumn', 1, 1, 150, 150)


# =============================================================================
# 5.2 - Create Columns
# =============================================================================

num_discretize = 5
ops.beamIntegration('Lobatto', 2, 2, 3)
ops.geomTransf('PDelta', 2)
opsdm.DiscretizeMember(1, 5, num_discretize, 'dispBeamColumn', 2, 2, 210, 210)
opsdm.DiscretizeMember(5, 9, num_discretize, 'dispBeamColumn', 2, 2, 250, 250)
opsdm.DiscretizeMember(9, 13, num_discretize, 'dispBeamColumn', 2, 2, 290, 290)
opsdm.DiscretizeMember(2, 6, num_discretize, 'dispBeamColumn', 2, 2, 220, 220)
opsdm.DiscretizeMember(6, 10, num_discretize, 'dispBeamColumn', 2, 2, 260, 260)
opsdm.DiscretizeMember(10, 14, num_discretize, 'dispBeamColumn', 2, 2, 2100, 2100)
opsdm.DiscretizeMember(3, 7, num_discretize, 'dispBeamColumn', 2, 2, 230, 230)
opsdm.DiscretizeMember(7, 11, num_discretize, 'dispBeamColumn', 2, 2, 270, 270)
opsdm.DiscretizeMember(11, 15, num_discretize, 'dispBeamColumn', 2, 2, 2110, 2110)
opsdm.DiscretizeMember(4, 8, num_discretize, 'dispBeamColumn', 2, 2, 240, 240)
opsdm.DiscretizeMember(8, 12, num_discretize, 'dispBeamColumn', 2, 2, 280, 280)
opsdm.DiscretizeMember(12, 16, num_discretize, 'dispBeamColumn', 2, 2, 2120, 2120)


# =============================================================================
# 6 - Element Tags
# =============================================================================

elements = ops.getEleTags()
beamTags = []
colTags = []

for element in elements:
    if element <= 200:
        beamTags.append(element)
    else:
        colTags.append(element)
    
        
# =============================================================================
# 7 - Plot the Model 
# =============================================================================
# opsv.plot_model(element_labels=1, node_labels=1)



# =============================================================================
# 8 - Apply loads for static analysis
# =============================================================================

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

g = 9.81
q_beam = 11330 
q_col = 1686.0937500000002
mass_endColumns= 3114.2250000000004 
mass_centralColumns = 5712.825000000002

#beam loads
for beamTag in beamTags:
    ops.eleLoad('-ele', beamTag, '-type', '-beamUniform', -q_beam, 0)

#columns loads
for colTag in colTags:
    ops.eleLoad('-ele', colTag, '-type', '-beamUniform', 0, -q_col)

# Concentrated masses at the nodes of the floor columns.
endColumns = [1, 4, 5, 8, 9, 12, 13, 16]
centralColuns = [2, 3, 6, 7, 10, 11, 14, 15]

for no in endColumns:
    ops.mass(no, mass_endColumns, mass_endColumns, 0.0)

for no in centralColuns:
    ops.mass(no, mass_centralColumns, mass_centralColumns, 0.0)
            
            
# =============================================================================
# 9 - Start of analysis generation
# =============================================================================

ops.system('BandGeneral')
ops.constraints('Transformation')
ops.numberer('RCM')
ops.test('NormDispIncr', 1.0e-8, 10, 3)
ops.algorithm('Newton')
ops.integrator('LoadControl', 0.1)
ops.analysis('Static')
ops.analyze(10)
print("Gravity Analysis Completed")


# =============================================================================
# 10 - Set the gravity loads to be constant & reset the time in the domain
# =============================================================================

ops.loadConst('-time', 0.0)


# =============================================================================
# 11 - Modal analys and set the Rayleigh Damping
# =============================================================================

eps_damping =  0.042
numEigen = 5 #Number of eigenvalores
eigenValues = ops.eigen('-genBandArpack',numEigen)  
w1 = np.sqrt(eigenValues[0])
w2 = np.sqrt(eigenValues[1])
w3 = np.sqrt(eigenValues[2])
T1_static = 2*3.1415/w1
T2_static = 2*3.1415/w2
T3_static = 2*3.1415/w3
alphaM = (1)*eps_damping*(2*w1*w3)/(w1*w3)
betaKcurr = (0)*2*eps_damping/(w1+w3)
betaKcomm = (1)*2*eps_damping/(w1+w3)
betaKinit = (0)*2*eps_damping/(w1+w3)
ops.rayleigh(alphaM, betaKcurr, betaKinit, betaKcomm)


#Reset the time in the domain
ops.loadConst('-time', 0.0)  #
   


# =============================================================================
# 12 - Preparing Earthquake Analysis
# =============================================================================

record = 'RSN1244_CHICHI_CHY101-E'
dt, nPts = ReadRecord.ReadRecord(record+'.at2', record+'.dat')


# =============================================================================
# 13 - Set an "absurdly" scaled timeseries
# =============================================================================

scale_factor = 2
ops.timeSeries('Path', 2, '-filePath', record+'.dat', '-dt', dt, '-factor', scale_factor*g, '-prependZero')
ops.pattern('UniformExcitation',  2,   1,  '-accel', 2)


# =============================================================================
# 14 - Just recorders
# =============================================================================

# Just record one node from each floor
# nodes_floors = [1,5,9,13]
# # Recorder: displacements of the floors
# ops.recorder('Node', '-file', '_Disp.out','-time', '-node', *nodes_floors, '-dof', 1,2,3, 'disp')

# # Recorder: relative acceleration of the floors
# ops.recorder('Node', '-file', '_AccelsREL.out','-time', '-node', *nodes_floors, '-dof', 1,2,3, 'accel')

# # Recorder: absolute acceleration of the floors
# ops.recorder('Node', '-file', '_AccelsABS_x.out', '-timeSeries', 2, '-node', *nodes_floors, '-dof', 1, 'accel')

# # Recorder: Obtain velocity of the floors
# ops.recorder('Node', '-file', '_Vel.out','-time', '-node', *nodes_floors, '-dof', 1,2,3, 'vel')



# =============================================================================
# 15 - Run Analysis
# =============================================================================

type_algorithm = {1:'NewtonLineSearch', 2: 'KrylovNewton' , 3:'ModifiedNewton' , 4: 'BFGS'}
aux_initial = ['KrylovNewton', 'ModifiedNewton']

ops.wipeAnalysis()
ops.system('BandGeneral')
ops.constraints('Transformation')
ops.numberer('RCM')

# set some variables
tFinal = nPts*dt
tCurrent = ops.getTime()
time = [tCurrent]

ok = 0
failure = 0

while tCurrent < tFinal and failure == 0:
    
    
    jj = 1
    ops.test('NormDispIncr', 1.0e-6, 1000, 2)
    ops.integrator('Newmark', 0.5, 0.25)
    ops.algorithm(type_algorithm[jj])
    ops.analysis('Transient')

    ok = ops.analyze(1, dt) #trocar para dt
    

    if ok != 0:
        while ok != 0 and jj < len(type_algorithm)+1:
            print('\n\n###############################################')
            print('FAILED !!', type_algorithm[jj], 'NormDispIncr', ' t= ', ops.getTime())
            print('###############################################\n\n')
            
            jj +=1
            
            if jj != len(type_algorithm)+1:
                
                if type_algorithm[jj] in aux_initial:
                    ops.algorithm(type_algorithm[jj], '-initial')
                    
                else:
                    ops.algorithm(type_algorithm[jj])
                    
                ok = ops.analyze(1, dt) 
    
    
    if ok != 0:
        jj = 1
        print('\n\n###############################################')
        print('    Try to reduce "dt" to "dt/10"')
        print('###############################################\n\n')
            
        while ok != 0 and jj < len(type_algorithm)+1:

            if type_algorithm[jj] in aux_initial:
                ops.algorithm(type_algorithm[jj], '-initial')
            else:
                ops.algorithm(type_algorithm[jj])
            
            
            ok = ops.analyze(1*10, dt/10)
            
            if ok != 0:
                print('\n\n###############################################')
                print('FAILED !!', type_algorithm[jj], dt/10, 'NormDispIncr', ' t= ', ops.getTime())
                print('###############################################\n\n')
                jj +=1
                
                
                
    if ok != 0:
        jj = 1
        print('\n\n###############################################')
        print('    Try to reduce "dt" to "dt/100"')
        print('###############################################\n\n')
            
        while ok != 0 and jj < len(type_algorithm)+1:

            if type_algorithm[jj] in aux_initial:
                ops.algorithm(type_algorithm[jj], '-initial')
            else:
                ops.algorithm(type_algorithm[jj])
            
    
            ok = ops.analyze(1*100, dt/100)
            
            if ok != 0:
                print('\n\n###############################################')
                print('FAILED !!', type_algorithm[jj], dt/100, 'NormDispIncr', ' t= ', ops.getTime())
                print('###############################################\n\n')
                jj +=1
    

    if ok == 0:
        tCurrent = ops.getTime()                
        time.append(tCurrent)

        print(f'{record[0:7]}_{scale_factor}', 'NormDispIncr', type_algorithm[jj], 't =',tCurrent, 'tF=', tFinal)
        
        
    else: 
        print('\n\n###############################################')
        print('ANALYSIS FAILED in t= ', ops.getTime())
        print('###############################################\n\n')
        failure = 1
        break
    
 
print(f"\nEnd of the Analysis: '{record[0:7]}'")
    

if ok == 0:
    feed_analysis = 'ANALYSIS SUCCESS'
else:
    feed_analysis = 'ANALYSIS FAILED'

print('\nAnalysis completed with:', feed_analysis)
    








