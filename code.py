"""
The following code is preamble:
"""

import mpmath as mpm
from mpmath import fdiv, exp, log, pi, power, sqrt

mpm.dps = 150

#def lnW(x,a,b):
#    return fdiv(3 - 8*a,3)*x - 2*log(b*x,2)

ranges = [
        (80, 90,     0.9850, 3.7065573*power(10,5), 5.216),
        (90, 100,    0.9850, 4.2572147*power(10,5), 4.831),
        (100, 110,   0.9850, 4.8890114*power(10,5), 4.513),
        (110, 120,   0.9850, 5.4574421*power(10,5), 4.264),
        (120, 130,   0.9850, 6.2949027*power(10,5), 4.032),
        (130, 140,   0.9850, 6.9404543*power(10,5), 3.855),
        (140, 150,   0.9850, 7.7137378*power(10,5), 3.696),
        (150, 160,   0.9850, 8.2591307*power(10,5), 3.572),
        (160, 170.2, 0.9850, 9.0996607*power(10,5), 3.448),
        (170.2, 500, 0.9850, 1.12*power(10,6),      3.337),
        (500, 1000,  0.9850, 6.23*power(10,6),      2.152),
        (1000, 1500, 0.9850, 2.12*power(10,7),      1.820),
        (1500, 2000, 0.9850, 6.47*power(10,7),      1.684),
        (2000, 2500, 0.9850, 1.73*power(10,8),      1.610),
        (2500, 3000, 0.9850, 2.56*power(10,8),      1.577),
        (3000, 481958, 0.9850, 5.76*power(10,8),    1.551),
        (481958, 6.7*power(10,12), 0.9850, 1.62*power(10,11), 1.448),
    ]

def BellottiTable(T):
    for T0, T1, alpha0, C, B in ranges:
        if T0 <= T < T1:
            return C, B
def CHL_verifier(k):
    return k*log(k*2.51949*power(10,11)*power(1 - fdiv(1,2.51949*power(10,11)),2))
def J(a,y):
    return fdiv(a*y,6) + log(a*y) + log(0.618)
def Rt(a,y):
    return fdiv(J(a,y) + 0.685 + 0.155*log(a*y),a*y*(0.04962 - fdiv(0.0196,J(a,y)+1.15)))
def Z_ve(x,m1):
    return min(Rt(m1,x)*log(m1*x),21.233,fdiv(53.989*power(log(x*m1),4/3),power(x*m1,1/3)))
def lnf(x,m1,option):
    if option == 1:
        C_Bellotti, B_Bellotti = BellottiTable(x*m1)
        return x*(1 - m1*B_Bellotti)
    if option == 2:
        return fdiv(3 - 8*m1,3)*x - 2*log(m1*x)
def upper_bound(x,m1,s1,k,option=2,auto=True):
    # This function returns the upper bound for the normalised difference of theta-functions
    # If this function outputs \delta < 1, then \theta(x + h) - \theta(x) > (1 - \delta) h..
    if auto == True:
        Z = Z_ve(x,m1)
    else:
        Z = auto
    h = k*exp(x*fdiv(k-1,k))
    delta = fdiv(k*exp(-fdiv(x,k)),x)
    if k >= 90:
        M = 6.391
        pwr = 0.9
    else:
        M = 1.26
        pwr = 0.2 
    alpha1 = 1 + 1.93378*power(10,-8)
    alpha2 = 1 + 1.936*power(10,-8)
    if option == 1:
        C_Bellotti, B_Bellotti = BellottiTable(x*m1)
        b1 = C_Bellotti
        b2 = 0
        b3 = 0
        b4 = 0
    if option == 2:
        b1 = 17.418
        b2 = 3
        b3 = 5.272
        b4 = 2
    t1 = fdiv(m1*x*exp(x*(s1 + m1 - 1)),pi) + 2*b3*power(m1*x,b4)*(power(m1*x,- fdiv(1,Z*m1)) - exp(x*(s1 - 1)))
    t2 = fdiv(2*b1*power(m1*x,b2)*x,lnf(x,m1,option))*(power(m1*x,-fdiv(lnf(x,m1,option),Z*m1*x)) - exp((s1-1)*lnf(x,m1,option)))
    t3 = fdiv(2*M,k)*fdiv((1 + k*exp(-fdiv(x,k)))*power(x + k*exp(-fdiv(x,k)),1-pwr) + power(x,1-pwr),exp(x*(m1 - fdiv(1,k))))
    t4 = (alpha1*sqrt(1 + k*exp(-fdiv(x,k))) - 0.999)*fdiv(1,k)*exp(x*(fdiv(1,k) - fdiv(1,2)))
    t5 = (alpha2*power(1 + k*exp(-fdiv(x,k)),fdiv(1,3)) - 0.885)*fdiv(exp(x*(fdiv(1,k) - fdiv(2,3))),k)
    return t1 + t2 + t3 + t4 + t5

"""
The following functions are used to calculate optimal choices for c and the computed ranges of x.
"""

def optimise(x,s1,k,option=2,auto=True):
    c = 1.00001
    if option == 1:
        while fdiv(x*c,k) < 80:
            c += 0.00001
    if option == 2:
        while fdiv(x*c,k) < log(3000175332800):
            c += 0.00001
    for q in [0.01,0.001,0.0001,0.00001]:
        while upper_bound(x,fdiv(c + q,k),s1,k,option,auto) < upper_bound(x,fdiv(c,k),s1,k,option,auto):
            c += q
    return c / k
def min_x(k,auto=True,refine=False):
    on = True
    mx = 100000
    while on == True:
        for o,s in [(2,0.675)]:
            a = optimise(mx,s,k,o,auto)
            ub = upper_bound(mx,a,s,k,o,auto)
        if ub < 1:
            on = False
        else:
            mx += 100000
    for q in [10000,1000,100,10,1]:
        on = True
        while on == True:
            if mx - q >= 1000:
                outcomes = []
                for o,s in [(2,0.675)]:
                    a = optimise(mx-q,s,k,o,auto)
                    outcomes.append((a*k,upper_bound(mx-q,a,s,k,o,auto),auto))
                #print(k, mx - q, outcomes)
                if outcomes[0][1] < 1:
                    mx -= q
                else:
                    on = False
            else:
                on = False
    if refine == True:
        for q in [10000,1000,100,10,1]:
            on = True
            while on == True:
                if mx - q >= 1000:
                    outcomes = []
                    for o,s in [(1,0.985)]:
                        a = optimise(mx-q,s,k,o,auto)
                        outcomes.append((a*k,upper_bound(mx-q,a,s,k,o,auto),auto))
                    #print(k, mx - q, outcomes)
                    if outcomes[0][1] < 1:
                        mx -= q
                        final_outcomes = outcomes
                    else:
                        on = False
                else:
                    on = False
    outcomes = []
    for o,s in [(1,0.985), (2,0.675)]:
        a = optimise(mx,s,k,o,auto)
        outcomes.append((a*k,upper_bound(mx,a,s,k,o,auto),o))
    for o in outcomes:
        if o[1] < 1:
            final_outcomes = o
    #print(mx, final_outcomes)
    return mx, final_outcomes

"""
The following code is used to compute the contents of Table 3.
"""

for kay in [100 - i for i in range(31)]:
    if kay >= 90:
        mx_kay = min_x(kay,True,False)
        print(kay, mx_kay[0], round(mx_kay[1][0],5))
    else:
        mx_kay = min_x(kay,True,True)
        mx_kay_False = min_x(kay,True,False)
        print(kay, mx_kay_False[0], round(mx_kay_False[1][0],5), mx_kay[0], round(mx_kay[1][0],5))


"""
The following code is used to compute the contents of Table 5 (i.e., the minimal zero-free regions).
"""

zed = 26.0
kay = 100
while kay >= 70:
    x0 = CHL_verifier(kay)
    t1, t2 = min_x(kay,zed)
    #print(zed, t1, t2)
    for q in [0.1,0.01,0.001]:
        while t1 > x0:
            zed -= q
            t1, t2 = min_x(kay,zed)
            #print(zed, t1, t2)
        zed += q
        t1, t2 = min_x(kay,zed)
    zed -= q
    t1, t2 = min_x(kay,zed)
    print(kay, round(zed,3), t1, round(t2[0],5))
    kay -= 1
