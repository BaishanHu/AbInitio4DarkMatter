import math
import DM_constants as dmc

def Fit_F(u, c):
    y = 0
    for i in range(len(c)):
        y += c[i] * u**i
    val = math.exp(-u/2.0)*y
    return val

# Fi[0]=F_Sigmapis; Fi[1]=F_Sigmapiv; Fi[2]=F_Sigmappis; Fi[3]=F_Sigmappiv;
def F2S(type, Fi, F_dSigmap, F_dSigmapp):
    if type=='p':
        data_temp = ( Fi[0]+(1.0+F_dSigmap)*Fi[1] )**2 + \
        ( Fi[2]+(1.0+F_dSigmapp)*Fi[3] )**2
    elif type=='n':
        data_temp = ( Fi[0]-(1.0+F_dSigmap)*Fi[1] )**2 + \
        ( Fi[2]-(1.0+F_dSigmapp)*Fi[3] )**2
    elif type=='00':
        data_temp = Fi[0]*Fi[0] + Fi[2]*Fi[2]
    elif type=='11':
        data_temp = ( (1.0+F_dSigmap)*Fi[1] )**2 + \
        ( (1.0+F_dSigmapp)*Fi[3] )**2
    elif type=='01':
        data_temp = (1.0+F_dSigmap)*Fi[1]*Fi[0] + \
        (1.0+F_dSigmapp)*Fi[3]*Fi[2]
        data_temp = data_temp * 2.0
    return data_temp

# Fi[0]=F_Sigmap_p; Fi[1]=F_Sigmap_n; Fi[2]=F_Sigmapp_p; Fi[3]=F_Sigmapp_n;
def F2S_pn(type, Fi, F_dSigmap, F_dSigmapp):
    if type=='p':
        data_temp = ( 2*Fi[0]+F_dSigmap*(Fi[0]-Fi[1]) )**2 + \
        ( 2*Fi[2]+F_dSigmapp*(Fi[2]-Fi[3]) )**2
    elif type=='n':
        data_temp = ( 2*Fi[1]-F_dSigmap*(Fi[0]-Fi[1]) )**2 + \
        ( 2*Fi[3]-F_dSigmapp*(Fi[2]-Fi[3]) )**2
    elif type=='00':
        data_temp = (Fi[0]+Fi[1])**2 + (Fi[2]*Fi[3])**2
    elif type=='11':
        data_temp = ( (1.0+F_dSigmap)*(Fi[0]-Fi[1]) )**2 + \
        ( (1.0+F_dSigmapp)*(Fi[2]-Fi[3]) )**2
    elif type=='01':
        data_temp = (1.0+F_dSigmap)*(Fi[0]+Fi[1])*(Fi[0]-Fi[1]) + \
        (1.0+F_dSigmapp)*(Fi[2]+Fi[3])*(Fi[2]-Fi[3])
        data_temp = data_temp * 2.0
    return data_temp

def dSigmap_1B_Klos(q):
    val = - 2.0*(q/dmc.LambdaA)**2
    return val

def dSigmapp_1B_Klos(q):
    val = -1.0*(dmc.gpipn*dmc.Fpi*q*q)/(dmc.M_NUCLEON*dmc.gA*(q*q+dmc.mpi*dmc.mpi))
    return val

def delta_2BC_Klos(rho, c3, c4, cD, q):
    # print('delta_2BC:', dmc.mpi, dmc.gA, dmc.LambdaChi, dmc.M_NUCLEON)
    rho = rho*dmc.HBARC**3; c3 = c3*1e-3; c4 = c4*1e-3
    kF = ( 3.0*math.pi*math.pi*rho/2.0 )**(1.0/3.0)

    if ( abs(q) > 1e-6 ):
        v_arccot = math.pi/2.0 - math.atan2( dmc.mpi**2 + q**2/4.0-kF**2, 2.0*dmc.mpi*kF )
        v_log = math.log( (dmc.mpi**2 + (kF-q/2.0)**2)/(dmc.mpi**2 + (kF+q/2.0)**2) )

        va1 = 8.0*kF*q*( 48.0*(kF**2+dmc.mpi**2)**2 + 32.0*(kF**2-3.0*dmc.mpi**2)*q**2-3.0*q**4)
        va2 = 768.0*dmc.mpi**3*q**3*v_arccot
        va3 = 3.0*( 16.0*(kF**2+dmc.mpi**2)**2 - 8.0*(kF**2-5.0*dmc.mpi**2)*q**2 + q**4 )
        va4 = 4.0*(kF**2+dmc.mpi**2) - q**2
        I1 = 1.0/(512.0*kF**3*q**3)*( va1+va2+va3*va4*v_log )

        vb1 = 8.0*kF*( 2.0*kF**2-3.0*dmc.mpi**2)*q
        vb2 = 24.0*dmc.mpi**3*q*v_arccot
        vb3 = 3.0*dmc.mpi**2*( 4.0*kF**2-q**2+4.0*dmc.mpi**2 )
        I2 = 1.0/(16.0*kF**3*q)*( vb1+vb2+vb3*v_log )

        vc1 = 8.0*kF*q*( 48.0*(kF**2+dmc.mpi**2)**2 - 32.0*kF**2*q**2-3.0*q**4)
        vc2 = 3.0*( 4.0*(kF**2+dmc.mpi**2) - q**2 )
        vc3 = 4.0*dmc.mpi**2 + (2.0*kF-q)**2
        vc4 = 4.0*dmc.mpi**2 + (2.0*kF+q)**2
        IP = -3.0/(512.0*kF**3*q**3)*( vc1+vc2*vc3*vc4*v_log )

        vd1 = 32.0*kF**3*q + 32.0*kF*dmc.mpi**2*q + 8.0*kF*q**3
        vd2 = 16.0*(kF**2+dmc.mpi**2)**2 + 8.0*(dmc.mpi**2-kF**2)*q**2 + q**4
        vd3 = math.log( (4.0*dmc.mpi**2 + (2.0*kF-q)**2)/(4.0*dmc.mpi**2 + (2.0*kF+q)**2) )
        Ic6 = -9.0/(128.0*kF**3*q)*( vd1+vd2*vd3 )
    else:
        v_arccot = math.atan2( kF, dmc.mpi )
        I1 = 1.0 - 3.0*(dmc.mpi/kF)**2 + 3.0*(dmc.mpi/kF)**3*v_arccot
        I2 = I1
        IP = 0.0
        Ic6 = 0.0

    da1 = 1.0/3.0*( c4+1.0/4.0/dmc.M_NUCLEON )*( 3.0*I2-I1 )
    da2 = 1.0/3.0*( -1.0*c3+1.0/4.0/dmc.M_NUCLEON )*I1
    da3 = (1.0+dmc.c6_hat)/12.0/dmc.M_NUCLEON*Ic6
    # da3 = (1.0+dmc.c6_hat)/4.0/dmc.M_NUCLEON*Ic6
    da4 = cD/( 4.0*dmc.gA*dmc.LambdaChi )
    delta_a1 = -1.0*rho/(dmc.Fpi**2)*( da1+da2-da3-da4 )

    db1 = -2.0*c3*q**2/(dmc.mpi**2+q**2)
    db2 = (c3+c4)/3.0*IP
    db3 = (1+dmc.c6_hat)/12.0/dmc.M_NUCLEON*Ic6
    # db3 = (1+dmc.c6_hat)/4.0/dmc.M_NUCLEON*Ic6
    delta_a1_P = rho/(dmc.Fpi**2)*( db1+db2-db3 )

    # if(abs(q) < 1e-6 ):
    #     print(kF, rho, c3, c4, I1, I2, delta_a1, delta_a1_P)

    if( (abs(c3) < 1e-6) and (abs(c4) < 1e-6) and (abs(cD) < 1e-6) ):
        delta_a1 = 0.0; delta_a1_P = 0.0

    return delta_a1, delta_a1_P

def dSigmap_1B(q):
    val = -1.0*q**2*dmc.r2A/dmc.HBARC/dmc.HBARC/6.0
    return val

def dSigmapp_1B(q):
    val = -1.0*(dmc.gpinn*dmc.Fpi*q*q)/(dmc.M_NUCLEON*dmc.gA*(q*q+dmc.mpi*dmc.mpi))
    return val

def delta_2BC(rho, c1, c3, c4, c6, cD, q):
    # print('delta_2BC:', dmc.mpi, dmc.gA, dmc.LambdaChi, dmc.M_NUCLEON)
    rho = rho*dmc.HBARC**3; c1 = c1*1e-3; c3 = c3*1e-3; c4 = c4*1e-3
    kF = ( 3.0*math.pi*math.pi*rho/2.0 )**(1.0/3.0)

    if ( abs(q) > 1e-6 ):
        v_arccot = math.pi/2.0 - math.atan2( dmc.mpi**2 + q**2/4.0-kF**2, 2.0*dmc.mpi*kF )
        v_log = math.log( (dmc.mpi**2 + (kF-q/2.0)**2)/(dmc.mpi**2 + (kF+q/2.0)**2) )

        va1 = 8.0*kF*q*( 48.0*(kF**2+dmc.mpi**2)**2 + 32.0*(kF**2-3.0*dmc.mpi**2)*q**2-3.0*q**4)
        va2 = 768.0*dmc.mpi**3*q**3*v_arccot
        va3 = 3.0*( 16.0*(kF**2+dmc.mpi**2)**2 - 8.0*(kF**2-5.0*dmc.mpi**2)*q**2 + q**4 )
        va4 = 4.0*(kF**2+dmc.mpi**2) - q**2
        I1 = 1.0/(512.0*kF**3*q**3)*( va1+va2+va3*va4*v_log )

        vb1 = 8.0*kF*( 2.0*kF**2-3.0*dmc.mpi**2)*q
        vb2 = 24.0*dmc.mpi**3*q*v_arccot
        vb3 = 3.0*dmc.mpi**2*( 4.0*kF**2-q**2+4.0*dmc.mpi**2 )
        I2 = 1.0/(16.0*kF**3*q)*( vb1+vb2+vb3*v_log )

        vc1 = 8.0*kF*q*( 48.0*(kF**2+dmc.mpi**2)**2 - 32.0*kF**2*q**2-3.0*q**4)
        vc2 = 3.0*( 4.0*(kF**2+dmc.mpi**2) - q**2 )
        vc3 = 4.0*dmc.mpi**2 + (2.0*kF-q)**2
        vc4 = 4.0*dmc.mpi**2 + (2.0*kF+q)**2
        IP = -3.0/(512.0*kF**3*q**3)*( vc1+vc2*vc3*vc4*v_log )

        vd1 = 32.0*kF**3*q + 32.0*kF*dmc.mpi**2*q + 8.0*kF*q**3
        vd2 = 16.0*(kF**2+dmc.mpi**2)**2 + 8.0*(dmc.mpi**2-kF**2)*q**2 + q**4
        vd3 = math.log( (4.0*dmc.mpi**2 + (2.0*kF-q)**2)/(4.0*dmc.mpi**2 + (2.0*kF+q)**2) )
        Ic6 = -9.0/(128.0*kF**3*q)*( vd1+vd2*vd3 )
    else:
        v_arccot = math.atan2( kF, dmc.mpi )
        I1 = 1.0 - 3.0*(dmc.mpi/kF)**2 + 3.0*(dmc.mpi/kF)**3*v_arccot
        I2 = I1
        IP = 0.0
        Ic6 = 0.0

    da1 = 1.0/3.0*( c4+1.0/4.0/dmc.M_NUCLEON )*( 3.0*I2-I1 )
    da2 = 1.0/3.0*( -1.0*c3+1.0/4.0/dmc.M_NUCLEON )*I1
    da3 = (1.0+c6)/12.0/dmc.M_NUCLEON*Ic6
    # da3 = (1.0+c6)/4.0/dmc.M_NUCLEON*Ic6
    da4 = cD/( 4.0*dmc.gA*dmc.LambdaChi )
    delta_a1 = -1.0*rho/(dmc.Fpi**2)*( da1+da2-da3-da4 )

    db1 = -2.0*(c3-2.0*c1)*(q*dmc.mpi)**2/(dmc.mpi**2+q**2)**2
    db2 = (c3+c4)/3.0*IP
    db3 = (1+c6)/12.0/dmc.M_NUCLEON*Ic6 - 2.0/3.0*c1*dmc.mpi**2/(dmc.mpi**2+q**2)*Ic6
    # db3 = (1+c6)/4.0/dmc.M_NUCLEON*Ic6 - 2.0/3.0*c1*dmc.mpi**2/(dmc.mpi**2+q**2)*Ic6
    db4 = q**2/(dmc.mpi**2+q**2)*( c3/3.0*(I1+IP)+(c4+1.0/4.0/dmc.M_NUCLEON)/3.0*(I1+IP-3.0*I2) )
    db5 = q**2/(dmc.mpi**2+q**2)*( cD/( 4.0*dmc.gA*dmc.LambdaChi ) )
    delta_a1_P = rho/(dmc.Fpi**2)*( db1+db2-db3-db4-db5 )

    # if(abs(q) < 1e-6 ):
    #     print(kF, rho, c3, c4, I1, I2, delta_a1, delta_a1_P)

    if( (abs(c1) < 1e-6) and (abs(c3) < 1e-6) and (abs(c4) < 1e-6) and (abs(cD) < 1e-6) ):
        delta_a1 = 0.0; delta_a1_P = 0.0

    return delta_a1, delta_a1_P
