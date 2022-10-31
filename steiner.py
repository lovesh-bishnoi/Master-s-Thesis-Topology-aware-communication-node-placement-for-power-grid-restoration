from math import *

#distance between two points
d = lambda x,y: ((x[0]-y[0])**2+(x[1]-y[1])**2)**0.5

#given the points, returns the sides 
s = lambda A,B,C : (d(B,C), d(C,A), d(A,B))

#given the sides, returns the angle
j = lambda a,b,c : acos((b*b+c*c-a*a)/(2*b*c))

#given the sides, returns secant of that angle
t = lambda a,b,c: 1/cos(j(a,b,c)-pi/6)

#given the sides and the Trilinear co-ordinates, returns the Cartesian co-ordinates
b = lambda A,B,C,p,q,r: [(p*A[i]+q*B[i]+r*C[i])/(p+q+r) for i in [0,1]] 

#this one checks if any of the angle is >= 2π/3 returns that point else computes the point
steiner_point = lambda A,B,C: A if j(*s(A,B,C)) >= 2*pi/3 else B if j(*s(B,C,A)) >= 2*pi/3 else C if j(*s(C,A,B)) >= 2*pi/3 else b(A,B,C,d(B,C)*t(*s(A,B,C)),d(C,A)*t(*s(B,C,A)),d(A,B)*t(*s(C,A,B)))

#  node 3,2,5
# print('{}'.format(steiner_point([550, 90], [150, 50], [300, 450])))
# print('{}'.format(steiner_point((550, 90), (150, 50), (300, 450) )))
# # [330.92512606555573, 184.59290966051807]

# #  node 3,1,5
# print('{}'.format(steiner_point([550, 90], [100,100], [300, 450])))
# # [304.94988731915925, 224.4005518349387]

# #  node 3,4,5
# print('{}'.format(steiner_point([550, 90], [70, 200], [300, 450])))
# [281.1576362100577, 290.73305842919336]