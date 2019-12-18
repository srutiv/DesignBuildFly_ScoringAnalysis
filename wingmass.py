"""
Script for estimating mass of a rectangular foam wing covered by balsa

Aron Zhao - az355@cornell.edu
November 6, 2019
"""


#gathering inputs
print('do not enter units')
b = float(input('enter wingspan (m): '))
thckB = float(input('enter balsa thickness (m): '))
dnsB = float(input('enter balsa density (kg/m^3): '))
dnsF = float(input('enter foam density (kg/m^3): '))
foilArea = float(input('enter airfoil area (m^2): '))
foilCircm = float(input('enter airfoil circumference (m): '))

#calculation for balsa
SA = foilCircm*b
Vb = SA*thckB
Mb = Vb*dnsB

#calculation for foam
Vf = foilArea*b
Mf = Vf*dnsF

#output
print(str(Mb + Mf) + ' kg')
input(press enter to exit)
