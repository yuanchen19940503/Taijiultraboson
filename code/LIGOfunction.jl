using PyPlot
# using QuadGK
using Cuba

# #"Heaviside Theta function"
# heaviside(t) = 0.5 * (sign(t) + 1)

using GWSC  
using DelimitedFiles:readdlm
using Interpolations
using SpecialFunctions
using Roots
using Memoize

function backup_matrix(array1 , fileName::String)

    if isdir("backup") == false
        mkdir("backup")
    end
    
    len = Int(sqrt(length(array1)))
    file = "backup/" * fileName

    open(file, "w") do io
        for i in 1:len
            for j in 1:len
                write(io, string(array1[i,j])*"  ")
            end
            write(io,"\n")
        end
    end
end

function backup3(array1, array2, array3 , fileName::String)

    if isdir("backup") == false
        mkdir("backup")
    end
    l1 = length(array1)
    l2 = length(array2)
    l3 = length(array3)
    if (l1 != l2) || (l1 != l3) || (l2 != l3)
        print("The two arrays do not have the same number of elements.
            I will quit now!")
        return
    end

    len = length(array1)
    file = "backup/" * fileName

    open(file, "w") do io
        for i in 1:len
            if array2[i] != Inf # throw the infinity numbers
                write(io, string(array1[i])*"  "*string(array2[i])*"  "*string(array2[i])*"\n")
            end
        end
    end
end

χs = range(0,1 , length=300);
H0_yr = 1/(14.018079396957893*1e9) # in the unit of yr , corresponding to H0=69.8 km/s/Mpc 
t0 = 13.8*1e9 #age of the Universe, yr
H0_SI = 2.26206394e-18 # SI unit

det_ligo = LIGO(name="LIGO_Design", TObs=4); #assuming a 4yr detection
det_ce = LIGO(name="CE", TObs=4); #assuming a 4yr detection
det_et = LIGO(name="ET", TObs=4); #assuming a 4yr detection
det_kagra = LIGO(name="KAGRA_Design", TObs=4);

function interLogxLogyIn(xs, ys; out=Inf)
    logxs = log10.(xs)
    logys = log10.(ys)
    
    logInter0 = interpolate((logxs,), logys, Gridded(Linear()))
    logInter = extrapolate(logInter0, out)
        
    # convert back to normal scale
    f -> 10^logInter(log10(f))
end
function interLogxIn(xs, ys; out=0)
    logxs = log10.(xs)
    
    logInter0 = interpolate((logxs,), ys, Gridded(Linear()))
    logInter = extrapolate(logInter0, out)
        
    # convert back to normal scale
    f -> logInter(log10(f))
end


γData = readdlm("backup/H1L1_orf.txt")
fsγ = γData[:, 1]
γs = γData[:, 2]
γLIGO = interLogxIn(fsγ, γs)
# lisa_points = readdlm("/home/yc/The-Stochastic-Gravitational-Wave-Background-from-the-non-linear-order-density-perturbations/delta_spectrum/code/backup/ΩPIslisa.txt");
# lisax = lisa_points[:,1];
# lisay = lisa_points[:,2];

# ligo_h = readdlm("/home/yc/BHB-LSD/code/backup/LIGO_design.txt");
# ASD_LIGO = interLogxLogyIn(ligo_h[:,1], ligo_h[:,2])
# Sn_LIGO(f) = ASD_LIGO(f)
# Ωn_LIGO(f) = (2*π^2/3/H0^2) * f^3 * Sn_LIGO(f)



# @memoize function LSD_Ωn(f)
#     if 1.0e4 <= f <= 3.0e5
#         2*π^2*f^3*LSD_Sn(f)/(3*H0^2)
#     end
# end
# BHMF_points = readdlm("/home/yc/Ultralight-DM/horizon/code/backup/BHMF.txt");
# BHMFx = BHMF_points[:,1];
# BHMFy = BHMF_points[:,2];

function ztot_int0(a)
    Ωm = 0.3
    ΩΛ = 0.7
    Ωr = 9.0e-5
    num = a
    den = (Ωm*a+ΩΛ*a^4.0 + Ωr)^(0.5)
    num / den 
end

function ztot(z)
    amin = 0 
    amax = (1.0+z)^(-1.0)
    z_tobeint(a) = (amax-amin)* ztot_int0(amin + (amax-amin)*a)
    function integrand(x, f)
        f[1] = z_tobeint(x[1])
    end
    result, err = cuhre(integrand ,  rtol=1e-1 , maxevals=Int(1e9))
    result[1]/H0_yr
end

zsgen = 10 .^ range(-5 , 7, length=30000)
ztot_inter = interLogxLogyIn(zsgen , ztot.(zsgen))

function ztot_inter0(z)
    if 0<= z <= 1.0e-5
        return 1.350811558678811e10
    elseif 1.0e-5 < z <1e7
        return ztot_inter(z)
    else
        return 0
    end
end



h_ns = readdlm("backup/hndata.txt");
h_nx = h_ns[:,1];
h_ny = h_ns[:,2];

h_ns_int = interLogxLogyIn(h_nx,h_ny)
function Sn_nsemo(f)
    if h_nx[1]<= f <= h_nx[end]
        Sn = h_ns_int(f) ^ (2.0) / f
        return Sn
    else 
        return 1
    end
end

Ωn_nsemo(f) = (2*π^2/3/H0_SI^2) * f^3 * Sn_nsemo(f)
ΩPInsemo = readdlm("backup/ΩPINSEMO.txt");
ΩPICE = readdlm("backup/ΩPICE.txt");
ΩPILSD = readdlm("backup/ΩPILSD.txt");
ΩPILIGO = readdlm("backup/ΩPILIGO.txt");
ΩPIET = readdlm("backup/ΩPIET.txt");

LSD_f = [1.0e4  , 2.0e4, 3.0e4  , 4.0e4  , 5.0e4,  6.0e4  , 1.0e5  , 3.0e5];
LSD_h = [2.9e-22, 7e-23, 4.0e-23, 2.5e-23, 2.0e-23,1.7e-23, 1.3e-23, 9e-24];
ASD_LSD = interLogxLogyIn(LSD_f,LSD_h)
Sn_LSD(f) = ASD_LSD(f)^2
Ωn_LSD(f) = (2*π^2/3/H0_SI^2) * f^3 * Sn_LSD(f)

function tL_D1(z) # tL in the unit of yr , this function is dt_L(z)/dz
    Δ = 0.3*(1.0+z)^3 +0.7
    1/H0_yr/sqrt(Δ)/(1+z)
end

function ztot_age(z) #approximated expression, see https://iopscience.iop.org/article/10.1088/1538-3873/aac1b2/pdf
    ΩΛ=0.7
    Ωm=0.3
    term1 = (ΩΛ/(1-ΩΛ))^0.5 * (1+z)^(-1.5)
    term2 = (ΩΛ/(1+z)^3/(1-ΩΛ)+1)^0.5
    1/H0_yr*2/3/sqrt(ΩΛ)*log(term1+term2)
end

function ttoz(t)  #change age t(yr) to redshift z
    ΩΛ=0.7
    Ωm=0.3
    fun = z -> ztot_age(z)-t
#     fun = z -> 1/H0_yr*2/3/sqrt(ΩΛ)*log((ΩΛ/(1-ΩΛ))^0.5 * (1+z)^(-1.5)+(ΩΛ/(1+z)^3/(1-ΩΛ)+1)^0.5) - t
    if t < 1e3
        return Inf
    elseif t <= ztot_age(0)
        return find_zero(fun , (0,1.0e5))
    else 
        return 0
    end
end

tsgen = 10 .^ range(-1 , 11, length=30000)
ttoz_inter = interLogxLogyIn(tsgen , ttoz.(tsgen))

function ttoz_inter0(t)
    if 1e-1<= t <= 1.0e11
        return ttoz_inter(t)
    elseif t<1e-1
        return Inf
    else
        return 0
    end
end

function tL(z) # tL in the unit of yr , lookback time analytical form, from https://doi.org/10.1088/1538-3873/aac1b2
    ΩΛ = 0.7
    2.0/(3.0*H0_yr*sqrt(ΩΛ))*log((1+ΩΛ^(-0.5))/((1.0+z)^(-3.0/2)+sqrt((1.0+z)^(-3.0)+(1-ΩΛ)/ΩΛ) ) )
end

@memoize function tLtoz(tL)  #change tL(yr) to redshift
    ΩΛ = 0.7    
    fun = z -> 2.0/(3.0*H0_yr*sqrt(ΩΛ))*log((1+ΩΛ^(-0.5))/((1.0+z)^(-3.0/2)+sqrt((1.0+z)^(-3.0)+(1-ΩΛ)/ΩΛ) ) ) - tL
    if tL*H0_yr < 0.9640
        return find_zero(fun , (0,1.0e5))
    else
        return 0.964099381639469  #this is tL(Inf)*H0_yr
    end
end

function ψ(z) #star formation rate in terms of redshift z  arXiv:1812.09622
    k = 0.178 # M⊙ /yr/Mpc³
    zm= 2.00
    a = 2.37
    b = 1.80
    n = k*a*exp(b*(z-zm))
    d = a - b +b*exp(a*(z-zm))
    n/d
end

function ϕ(Mstar) #initial mass function =a*M^(-2.35), a is the normalized coefficient
    a = 0.171636  #normalized in the mass range [0.1-100]M⊙ , ∫ ϕ(M) M dM = 1
    if 0.1<Mstar<100
        return a * Mstar^(-2.35)
    else
        return 0
    end
end

function Z_int0(z)
    ΩM=0.3
    Ωλ=0.7
    Ez = sqrt(ΩM*(1.0+z*1.0)^3.0+Ωλ)
    sfr = 0.015*(1.0+1.0*z)^(2.7) / (1.0 + ((1.0+1.0*z)/2.9)^(5.6) )
    97.8e10*sfr/69.8/Ez/(1+z)
end

function Z_int1(z)
    zmin = z
    zmax = 20
    Z_tobeint(z) = (zmax-zmin)* Z_int0(zmin + (zmax-zmin)*z)
    function integrand(x, f)
        f[1] = Z_tobeint(x[1])
    end
    result, err = cuhre(integrand ,  rtol=1e-3 , maxevals=Int(1e9))
    result[1]
end

function Z(z)#   arXiv: 1602.04531
    y = 0.019
    R = 0.27
    Ωb=0.045
    h0 = 0.7
    ρb = 2.77*1e11*Ωb*h0^2
    if z>20 || z<0
        return 0
    else
        logz = 0.5 + log10(y*(1-R)/ρb*Z_int1(z))
        return 10^(1.0*logz)
    end
end

@memoize function g_inv(mbh, z) #  arXiv: 1110.1726
    Zsun = 0.0196 #arXiv: 1703.10834
    if 3<=mbh<12
        fun = Mstar -> 1.1+0.2*exp((Mstar-11.0)/4.0) - (2.0 + Z(z)/Zsun)*exp(0.4*(Mstar-26.0)) - mbh
        find_zero(fun , (0,30))
    elseif 12<=mbh<=50
#         fun1= Mstar -> min(33.35 +(4.75+1.25*Z(z)/0.02)*(Mstar-34) , Mstar-sqrt(Z(z)/0.02)*(1.3*Mstar-18.35))-mbh
#         find_zero(fun1, (0 , 500))
        Mmax1 = ( -33.35+1.0*mbh+34.0*(4.75+1.25*Z(z)/Zsun) ) / (4.75+1.25*Z(z)/Zsun)
        Mmax2 = ( 1.0*mbh-18.35*sqrt(Z(z)/Zsun) ) / (1.0-1.3*sqrt(Z(z)/Zsun))
        return max(Mmax1, Mmax2)
    else
        return 0
    end
end

function lifetime_yr(M) #arXiv: astro-ph/0110697
    a0 = 9.785
    a1 = -3.759
    a2 = 1.413
    a3 = -0.186
    x = log10(M)
    10.0^(a0+a1*x+a2*x^2+a3*x^3)
end
    

function ndot_ligo(M , z)
    ms = g_inv(M , z) #mass of the star
    lifetime_z  = tLtoz(lifetime_yr(ms))
    if z<=lifetime_z
        return 0
    else
        return ψ(z-lifetime_z)*ϕ(ms)
    end
end

ndot_ligo_Ms = range(3,50,length=500)
ndot_ligo_zs = range(0.001,20,length=500)

# @time ndot_ligos = [ndot_ligo(M,z) for M in ndot_ligo_Ms, z in ndot_ligo_zs] #only need to run once
# backup_matrix(ndot_ligos,"ndot_ligos_itpuse.txt")

ndot_ligos = readdlm("backup/ndot_ligos_itpuse.txt");
ndot_ligo_itp = interpolate(ndot_ligos, BSpline(Cubic(Line(OnGrid()))))

ndot_ligo_sitp = scale(ndot_ligo_itp, ndot_ligo_Ms, ndot_ligo_zs)

ndot_ligo_inter = extrapolate(ndot_ligo_sitp,0);  #  1000 times faster!!!


function ndot_PTA(M , z)
    if 0<z<3
        return 0.005*(M/3.0e6)^(-0.3)/M/log(10.0)/t0  #1/Mpc^{-3}/yr/Msun
    else 
        return 0
    end
end

function ndot_PBH(fpbh , M , z)
    Mmin = 3.0
    Mmax = 50.0
    if 0<z<10
        if Mmin <= M <= Mmax
            massfunc = 1.0/2 *(Mmin^(-1/2) - Mmax^(-1/2))^(-1.0)* M^(-3/2.0)
            ρDM = 2.6e16 # Mdot/Mpc^3
            dndM = fpbh * ρDM * massfunc /M
            ndot = dndM / t0
            return ndot
        else
            return 0
        end
    else 
        return 0
    end
end





function μ(ms)
    ms*(7.48481*1.0e9)
end

function fs(ms)
    μ = ms*(7.48481*1.0e9)
    μ/π/(1.56235e-13)/365/86400
end


function isgrow(χ, ms , M ,m=1)
    Mfinal = Mf(χ, ms , M , m)
    MSmax = M - Mfinal
    l=m
    Γ = Γnlm(M , χ, ms , l+1 , l , m) #usually n = l+1
    tauinst = 1/Γ
    taugw = τGW(M , χ, ms  , l , m)
    if MSmax<=0 || Mfinal<=0 || Γ<=0 || (1/Γ>t0) #|| (taugw<=tauinst)
        return 0
    else
        return 1
    end
end

function Mf(χ, ms , M,m=1)
    μ = ms*(7.48481*1.0e9)
    x = μ*M
    if  m^6 - 16*m^2*x^2*(m-x*χ)^2  < 0
        return 0
    end
    Mf = M*(m^3 - sqrt( m^6 - 16*m^2*x^2*(m-x*χ)^2) ) / (8*x^2*(m-χ*x))
    if imag(Mf)!=0 || Mf<0
        return 0
    else
        return Mf
    end
end

@memoize function χcrit(M  , ms , m)
    isgrows = isgrow.(χs , ms , M , m)
    result = 0
    if maximum(isgrows)>0
        result = minimum(χs[map(x-> x>0 ,isgrows)])
    end
    if result>0
        return result
    else
        return NaN
    end
end

function Cnl(n,l) #checked
    term1 = 2.0^(4*l+1.0)
    term2 = n^(2.0*l+4.0)
    term0 = factorial(n+l)/factorial(n-l-1.0)
    term3 = factorial(l)
    term4 = factorial(2.0*l)*factorial(2.0*l+1)
    result = term1 / term2 *term0 *(term3/term4)^2.0
    return result
end

function dEdtCnl(n,l)
    num = 16.0^(l+1.0)*l*(2*l-1)*gamma(2*l-1.0)^2*gamma(l+n+1)^2.0
    den = n^(4*l+8.0)*(l+1.0)*gamma(l+1.0)^4.0 * gamma(4.0*l+3.0)*gamma(n-l)^2.0
    num/den
end


function dEdt(M,ms,n,l)
    μ = ms*(7.48481*1.0e9) # in the unit of 1/M⊙
    dEdtCnl(n,l)*(μ*M)^(4*l+10.0)
end

# dEdt(5,1e-12 , 2 ,1)

function τGW(M , χ , ms , l , m)
    μ = ms*(7.48481*1.0e9)
    x = M*μ
    det = m^6 - 16*m^2*x^2*(m-x*χ)^2
    if det < 0
        return 0
    end
    Mf = M*(m^3 - sqrt( m^6 - 16*m^2*x^2*(m-x*χ)^2) ) / (8*x^2*(m-χ*x))
    MSmax = M-Mf
    if Mf<=0 || MSmax<=0
        return 0
    else
        result = Mf*(dEdt(M,ms,l+1,l)*MSmax/Mf)^(-1.0) * 1.57092e-13 # change the unit of M⊙ to yr
    end
end

function τGW_appro(M,ms,χ)
    5.0e11 *((M/1.0e6)^14.0*(ms/1.0e-17)^15.0*χ)^(-1.0)
end
# τGW_appro(5,1e-12,0.2),τGW(5,0.2,1e-12,1,1)

function glm(M,χ,ms,l,m)
    a     = χ * M
    μ = ms*(7.48481*1.0e9) #unit = 1/M⊙
    rplus = M + sqrt(M^2-a^2) #unit = M⊙
    α = μ*M
    n = l+1 #dominant mode n>=l+1
    fnl = - 6.0/(2*l+1)+2/n
    hl  = 16.0/(2*l*(2*l+1)*(2*l+2))
    ν   = n + fnl*α^2 + hl * α^3.0
    if α<=ν
        ω = μ*sqrt(1-α^2/ν^2)
    else
        return 0
    end
    result = 1
    for k = 1:l
        result=result * (k^2.0 * (1.0-χ^2.0) + (χ*m - 2.0*rplus * ω)^2.0) 
    end
    return result
end

# function Γnlm(M , χ, ms , n,l,m) #usually n = l+1
#     μ = ms*(7.48481*1.0e9) # in the unit of 1/M⊙
#     a     = χ * M
#     rplus_tilde = (M + sqrt(M^2-a^2))/M
#     rplus = (M + sqrt(M^2-a^2))
#     ΩH = a/2/M/rplus  # in the unit of 1/M⊙
#     α = μ*M
#     α_SI = α*2.11456#  *1.0e15  in the unit of kg⁻²
#     G_SI = 6.67408 # * 1.0e-11
#     c_SI = 2.99792458  # *1.0e8
#     hbar_SI = 1.0545718  # * 1.0e-34
#     result = 2*rplus_tilde*Cnl(n,l)*glm(M,χ,ms,l,m)*(m*ΩH-μ)*α_SI^(4*l+5) / 1.9891e30 #change M⊙ to kg
#     coe_unittoyr = G_SI^(-2*(3+2*l))*hbar_SI^(5+4*l)*c_SI^(4*(2+l)) * 31536000.0
#     #coe_unittoyr is in the unit of yr⁻¹
#     float_num = 15.0*(4*l+5)+(-11.0)*-2*(3+2*l)+(-34.0)*(5+4*l)+8.0*(4*(2+l))
# #     float_num is α_SI^(4*l+5)*G_SI^(-2*(3+2*l))*hbar_SI^(5+4*l)*c_SI^(4*(2+l))
#     result = result * coe_unittoyr *10.0^(float_num)
# end

function Γnlm(M, χ , ms , n , l , m)
    μ = ms*(7.48481*1.0e9)
    a     = χ * M
    rplus = (M + sqrt(M^2-a^2))    
    α = μ*M
    α_SI = α*2.11456#  *1.0e15  in the unit of kg⁻²
    G_SI = 6.67408 # * 1.0e-11
    c_SI = 2.99792458  # *1.0e8
    hbar_SI = 1.0545718  # * 1.0e-34
    ΩK = a/(a^2.0+rplus^2.0)
    p = n-l-1
    κK = (rplus^2.0-a^2.0)/(2.0*rplus*(rplus^2.0+a^2.0))
    term1 = (1.0*l)^(-4.0*l-9/2.0+p*1.0)
    term2 = 2.0^(2.0*l+1.0-1.0*p)*sqrt(π)*factorial(p)
    term3 = α^(4.0*l+5.0)
    termsinh = sinh(π*(l*ΩK-μ)/κK)
    termexp = -2/κK*(l*ΩK-α/rplus)*atan(ΩK/κK)-2.0*(1.0-l+p)
    result = term1/term2*term3*termsinh*exp(termexp)*M / 1.9891e30 #change M⊙ to kg
    coe_unittoyr = G_SI^(-2.0*(3.0+2.0*l))*hbar_SI^(5.0+4.0*l)*c_SI^(4.0*(2.0+l)) * 31536000.0
    float_num = 15.0*(4*l+5)+(-11.0)*-2*(3+2*l)+(-34.0)*(5+4*l)+8.0*(4*(2+l))
    result = result * coe_unittoyr *10.0^(float_num)
end

@memoize function EGW0(M , χ , ms,  tf , m=1) #only consider m=1
    μ = ms*(7.48481*1.0e9)
    ωR = μ    
    χc = χcrit(M  , ms , m)
    if χ >= χc
        Mfinal = Mf(χ, ms , M , m)
        MSmax = M - Mfinal
        l = m
        Γ = Γnlm(M , χ, ms ,l+1.0,l,m)
        if MSmax<=0 || Mfinal<=0 || Γ<=0 || 1/Γ>=t0
            return 0
        else
            τ = τGW(M , χ , ms ,l, m)
            if τ<0
                return 0
            else
#                 Δt = min(τ , t0)
                Δt = ztot_inter0(0) - tf
                return  MSmax*Δt/(τ+Δt)
            end
        end
    else
        return 0
    end
end

@memoize function EGW(M , χ ,  ms , tf ) #up to (n,l,m)=(2,1,1),(3,2,2)...(m+1,m,m),maximum m =10.0
    μ = ms*(7.48481*1.0e9)
    ωR = μ
    Egw = 0
    χrest = χ
    Mrest = M
    Δt = ztot_inter0(0) - tf
    i=1.0
    for i=1.0:10.0
        if χrest > χcrit(M  , ms , i)
            if (0 < 1/Γnlm(Mrest , χrest, ms , i+1 , i , i) < Δt) && Mrest>0 && χrest>0
                Egw = Egw + EGW0(Mrest , χrest , ms ,  tf , i)
                χrest = χcrit(Mrest  , ms , i)
                Mrest = Mf(χrest , ms , Mrest , i)
            end
        end
    end
    return Egw
end

function Ωgw_int0(f, M  , χ , ms , model) #the input f is in the unit of Hz 
    f = f*365*86400.0
    μ = ms*(7.48481*1.0e9)  #in the unit of 1/M⊙
    fs = μ/π/(1.56235e-13) #in the unit of 1/yr
    zmax = 20
    fmin = fs/(zmax+1.0)
    if fmin < f < fs
        z = fs/f - 1.0
        tf = ztot_inter0(z)
        if model == "old"
            Egw = EGW0( M , χ , ms, tf)
        elseif model == "new"
            Egw = EGW(M , χ , ms, tf)
        end
        ρc = 1.35174*10^11  # 3H0^2/(8πG)
#         ndot = ndot_ligo_inter(M , z)
#         ndot = ndot_PBH(M,z)
        
    ndot = ndot_PTA(M , z)
    return 1/ρc*tL_D1(z)*ndot*Egw 
    else
        return 0
    end
end

function Ωgw_int1(f , M , χ , ms , χmin , χmax, model )
#     Mmin = 1.0e8
#     Mmax = 1.0e15
    Mmin = 1e4
    Mmax = 1e7
    (Mmax-Mmin) * Ωgw_int0(f, Mmin + (Mmax-Mmin)*M  , χmin + (χmax-χmin)*χ , ms  , model)
end

@memoize function Ωgw(f , ms , χmin , χmax , model)
    function integrand(x, func)
        func[1] = Ωgw_int1(f , x[2]  , x[1] , ms , χmin , χmax , model)
    end
    result, err = cuhre(integrand , 2 , 1 ,  rtol=1e-1 , maxevals=Int(1e9))
    result[1]
end

zcutoff = 100.0

function fres(ms)
    fmin = log10(fs(ms)/(1+zcutoff))
    fmax = log10(fs(ms)*1.0)
    10 .^ range(fmin , fmax , length =50)
end

function SNR_int0(f, ms , χmin , χmax , det , TObs = 4.0)
    
    if det == "LSD"
        Pn = Sn_LSD(f) / 5.0
    elseif det =="NSEMO"
        Pn = Sn_nsemo(f) / 5.0
    else
        detector0 = LIGO(name=det, TObs=TObs);
        Pn = detector0.Sn(f) / 5.0
    end
    Ωs = Ωgw.(fres(ms) , ms , χmin , χmax , "new") + ΩgwSBH.(fres(ms),ms,"new")
#     Ωs = ΩgwSBH.(fres(ms),ms,"new")
    Ωhs = interpolate((fres(ms),), Ωs, Gridded(Linear()))
    if 1.0*fres(ms)[1] < 1.0*f < 1.0*fres(ms)[end]
        if det == "LIGO"
            Rf = γLIGO(f)/5  #R_eff(f) = Gamma(f)
        else
            Rf = 1.0/5   #Gaamma(f) = 1
        end
        Ωh = Ωhs(f)
        Sh = 3.0*H0_SI^2.0 * Ωh / 2.0 /f^3.0 / π^2.0
        num = Sh^2.0 * Rf^2.0
        den = (1.0/25+Rf^2.0)*Sh^2.0+Pn^2.0 + 2.0*Sh*Pn/5.0
        Tobs = TObs*365.0*86400.0
        return num/den*Tobs
    else
        return 0
    end
end

@memoize function SNR_ms(ms  ,χmin , χmax , det , TObs = 4.0)
    if det == "LSD"
        fmin = 1.0e4
        fmax = 3.0e5
    elseif det =="NSEMO"
        fmin = h_nx[1]
        fmax = h_nx[end]
    else
        detector = LIGO(name=det, TObs=TObs);
        fmin = detector.fPlotRange[1]
        fmax = detector.fPlotRange[2]
    end
    
    SNR_ms_tobeint(f,ms  ,χmin , χmax, det , TObs) = (fmax-fmin)* SNR_int0(fmin + (fmax-fmin)*f , ms  ,χmin , χmax , det , TObs)
    
    function integrand(x, f)
        f[1] = SNR_ms_tobeint(x[1],  ms  ,χmin , χmax, det , TObs)
    end
    result, err = cuhre(integrand ,  rtol=1e-1 , maxevals=Int(1e9))
    result[1]^(0.5)
end

function Rdenpbh(z , M1, M2 , α, Mpbh_min , Mpbh_max, fpbh) #in the unit of Gpc^-3 yr^-1, merger rate for PBHs
    Mmin = Mpbh_min
    Mmax = Mpbh_max
    t = ztot_inter0(z)
    p(M) = (1.0+α)/(Mmax^(1.0+α) -Mmin^(1.0+α))* M^(1.0*α)  # power law
#     σ=0.5
#     mc=3
#     p(M) = 1/sqrt(2*pi)/σ/M*exp(-log(M/mc)^2.0 / 2 / σ^2)#lognormal
    if (Mmin<=M1<=Mmax) && (Mmin<=M2<=Mmax)
        if M1==M2
            p1 = p(M1)/2
            p2 = p1
        else
            p1 = p(M1)
            p2 = p(M2)
        end
    else
        return 0
    end
    f = fpbh*0.85
    σeq = 0.005
    term1 = 0.0039*(t/t0)^(-34.0/37)*f^2*(f^2+σeq^2)^(-21.0/74) 
    term2 = min(p1/M1,p2/M2)*(p1/M1+p2/M2)
    term3 = (M1*M2)^(3.0/37) *(M1+M2)^(36.0/37)
    return term1*term2*term3
end

function ΩgwPBH_int0(fpbh , f, M1 , M2  ,α , ms , Mpbh_min , Mpbh_max, model) #the input f is in the unit of Hz 
    #the ogw from PBH remnants
#     χ = 1.0 - 1.0/(1.0+exp(M-30.0))
    if M1>M2
        return 0
    else
        ν=M1*M2/(M1+M2)^2.0
        M = M1+M2-(M1+M2)*((1-sqrt(8/9))*ν - 4*ν^2*(0.19308+sqrt(8/9)-1))
        χ = ν*(2*sqrt(3)-3.5171*ν+2.5763*ν^2.0)
#         χ = 0.7
#         M = M1 + M2 - 5.7e-2*(M1*M2)/(M1+M2) #Eq.27 in arXiv:1812.09622 the mass of remnant BHs
        f = f*365*86400.0
        μ = ms*(7.48481*1.0e9)  #in the unit of 1/M⊙
        fs = μ/π/(1.56235e-13) #in the unit of 1/yr
        fmin = fs/(zcutoff+1.0)
        if fmin < f < fs
            z = fs/f - 1.0
            tf = ztot_inter0(z)
            if model == "old"
                Egw = EGW0( M , χ , ms , tf)
            elseif model == "new"
                Egw = EGW(M , χ , ms , tf)
            end
            ρc = 1.35174*10^11  # 3H0^2/(8πG)
#             ndot = ndot_PBH(fpbh, M,z)
            Mmin = Mpbh_min
            Mmax = Mpbh_max

            Rate_den = Rdenpbh(z , M1, M2 , α , Mpbh_min , Mpbh_max, fpbh)
            return 1/ρc*tL_D1(z)*Rate_den*Egw
#         return 1/ρc*tL_D1(z)*ndot*Egw
        else
            return 0
        end
    end
end

function ΩgwPBH_int1(fpbh , f , M1 , M2 , α,  ms , Mpbh_min , Mpbh_max , model )
    Mmin = 1.0 #1.9715
    Mmax = 5.0 #9.8575
    (Mmax-Mmin)^(2.0) * ΩgwPBH_int0(fpbh , f, Mmin + (Mmax-Mmin)*M1, Mmin + (Mmax-Mmin)*M2  , α ,  ms ,  Mpbh_min , Mpbh_max  , model)
end

@memoize function ΩgwPBH(fpbh , f , ms, Mpbh_min , Mpbh_max, α , model)
    function integrand(x, func)
        func[1] = ΩgwPBH_int1(fpbh , f , x[1] , x[2] , α,  ms ,Mpbh_min , Mpbh_max,  model)
    end
    result, err = cuhre(integrand , 2 , 1 ,  rtol=1e-1,  maxevals=Int(1e9) )
    result[1]
end

###############################    log-normal    #######################################
function Rdenpbhlog(z , M1, M2 , mc, σ , fpbh) #in the unit of Gpc^-3 yr^-1, merger rate for PBHs
    t = ztot_inter0(z)
    p(M) = 1/sqrt(2*pi)/σ/M*exp(-log(M/mc)^2.0 / 2 / σ^2)#lognormal
    if M1==M2
        p1 = p(M1)/2
        p2 = p1
    else
        p1 = p(M1)
        p2 = p(M2)
    end
    f = fpbh*0.85
    σeq = 0.005
    term1 = 0.0039*(t/t0)^(-34.0/37)*f^2*(f^2+σeq^2)^(-21.0/74) 
    term2 = min(p1/M1,p2/M2)*(p1/M1+p2/M2)
    term3 = (M1*M2)^(3.0/37) *(M1+M2)^(36.0/37)
    return term1*term2*term3
end

function ΩgwPBHlog_int0(fpbh , f, M1 , M2  , ms , mc, σ, model) #the input f is in the unit of Hz 
    #the ogw from PBH remnants
#     χ = 1.0 - 1.0/(1.0+exp(M-30.0))
    if M1>M2
        return 0
    else
        ν=M1*M2/(M1+M2)^2.0
        M = M1+M2-(M1+M2)*((1-sqrt(8/9))*ν - 4*ν^2*(0.19308+sqrt(8/9)-1))
        χ = ν*(2*sqrt(3)-3.5171*ν+2.5763*ν^2.0)
#         χ = 0.7
#         M = M1 + M2 - 5.7e-2*(M1*M2)/(M1+M2) #Eq.27 in arXiv:1812.09622 the mass of remnant BHs
        f = f*365*86400.0
        μ = ms*(7.48481*1.0e9)  #in the unit of 1/M⊙
        fs = μ/π/(1.56235e-13) #in the unit of 1/yr
        fmin = fs/(zcutoff+1.0)
        if fmin < f < fs
            z = fs/f - 1.0
            tf = ztot_inter0(z)
            if model == "old"
                Egw = EGW0( M , χ , ms , tf)
            elseif model == "new"
                Egw = EGW(M , χ , ms , tf)
            end
            ρc = 1.35174*10^11  # 3H0^2/(8πG)
#             ndot = ndot_PBH(fpbh, M,z)
            Rate_den = Rdenpbhlog(z , M1, M2 , mc, σ , fpbh)
            return 1/ρc*tL_D1(z)*Rate_den*Egw
#         return 1/ρc*tL_D1(z)*ndot*Egw
        else
            return 0
        end
    end
end

function ΩgwPBHlog_int1(fpbh , f, M1 , M2  , ms , mc, σ, model)
    Mmin = 0.3
    Mmax = 10.0 
    (Mmax-Mmin)^(2.0) * ΩgwPBHlog_int0(fpbh , f, Mmin + (Mmax-Mmin)*M1, Mmin + (Mmax-Mmin)*M2  ,  ms , mc, σ, model)
end

@memoize function ΩgwPBHlog(fpbh , f , ms , mc, σ, model)
    function integrand(x, func)
        func[1] = ΩgwPBHlog_int1(fpbh , f , x[1] , x[2] ,  ms , mc, σ, model)
    end
    result, err = cuhre(integrand , 2 , 1 ,  rtol=1e-1,  maxevals=Int(1e9) )
    result[1]
end
########################################################################################

#####################################################################################################################
function ΩgwSBH_int0(f, M1 , M2  , ms , model) #the input f is in the unit of Hz 
    #the ogw from SBH remnants
    if 3<=M2<=M1 && M1+M2<=100
#         χ = 0.7
        ν=M1*M2/(M1+M2)^2.0
        M = M1+M2-(M1+M2)*((1-sqrt(8/9))*ν - 4*ν^2*(0.19308+sqrt(8/9)-1))
        χ = ν*(2*sqrt(3)-3.5171*ν+2.5763*ν^2.0)
#         M = M1 + M2 - 5.7e-2*(M1*M2)/(M1+M2) #Eq.27 in arXiv:1812.09622 the mass of remnant BHs
        f = f*365*86400.0
        μ = ms*(7.48481*1.0e9)  #in the unit of 1/M⊙
        fs = μ/π/(1.56235e-13) #in the unit of 1/yr
        fmin = fs/(20+1.0)
        if fmin < f < fs
            z = fs/f - 1.0
            tf = ztot_inter0(z)
            if model == "old"
                Egw = EGW0( M , χ , ms, tf)
            elseif model == "new"
                Egw = EGW(M , χ , ms , tf)
            end
            ρc = 1.35174*10^11  # 3H0^2/(8πG)
            Rate_den = Rbirth(z , M1, M2)
            return 1/ρc*tL_D1(z)*Rate_den*Egw
        else
            return 0
        end
    else
        return 0
    end
end

function ΩgwSBH_int1(f , M1 , M2 ,   ms , model )
    Mmin = 3.0  #5.9145
    Mmax = 50.0 #98.575
    (Mmax-Mmin)^(2.0) * ΩgwSBH_int0(f, Mmin + (Mmax-Mmin)*M1, Mmin + (Mmax-Mmin)*M2  ,  ms ,  model)
end

@memoize function ΩgwSBH(f , ms, model)
    function integrand(x, func)
        func[1] = ΩgwSBH_int1( f , x[1] , x[2] , ms , model)
    end
    result, err = cuhre(integrand , 2 , 1 ,  rtol=1e-1,  maxevals=Int(1e9))
    result[1]
end

#####################################################################################################################
function SNRPBH_int0(f, ms , Mpbh_min , Mpbh_max, α , det , TObs = 4.0)
    
    if det == "LSD"
        Pn = Sn_LSD(f) / 5.0
    elseif det =="NSEMO"
        Pn = Sn_nsemo(f) / 5.0
    else
        detector0 = LIGO(name=det, TObs=TObs);
        Pn = detector0.Sn(f) / 5.0
    end
    
    Ωhs = interpolate((fres(ms),), ΩgwPBH.(1e-3 , fres(ms), ms, Mpbh_min , Mpbh_max, α , "old"), Gridded(Linear()))
    if 1.0*fres(ms)[1] < 1.0*f < 1.0*fres(ms)[end]
        if det == "LIGO"
            Rf = γLIGO(f)/5  #R_eff(f) = Gamma(f)
        else
            Rf = 1.0/5   #Gaamma(f) = 1
        end
        Ωh = Ωhs(f)
        Sh = 3.0*H0_SI^2.0 * Ωh / 2.0 /f^3.0 / π^2.0
        num = Sh^2.0 * Rf^2.0
        den = (1.0/25+Rf^2.0)*Sh^2.0+Pn^2.0 + 2.0*Sh*Pn/5.0
        Tobs = TObs*365.0*86400.0
        return num/den*Tobs
    else
        return 0
    end
end

@memoize function SNRPBH_ms(ms  , Mpbh_min , Mpbh_max, α , det , TObs = 4.0)
    if det == "LSD"
        fmin = 1.0e4
        fmax = 3.0e5
    elseif det =="NSEMO"
        fmin = h_nx[1]
        fmax = h_nx[end]
    else
        detector = LIGO(name=det, TObs=TObs);
        fmin = detector.fPlotRange[1]
        fmax = detector.fPlotRange[2]
    end
    
    SNRPBH_ms_tobeint(f, ms  ,Mpbh_min , Mpbh_max, α , det , TObs) = (fmax-fmin)* SNRPBH_int0(fmin + (fmax-fmin)*f , ms  ,Mpbh_min, Mpbh_max, α ,det , TObs)
    
    function integrand(x, f)
        f[1] = SNRPBH_ms_tobeint(x[1],  ms  ,Mpbh_min , Mpbh_max, α , det , TObs)
    end
    result, err = cuhre(integrand ,  rtol=1e-2 , maxevals=Int(1e9))
    result[1]^(0.5)
end

######################################################################################################################################
function SNRPBHlog_int0(f, ms ,  mc, σ , det , TObs = 4.0)
    if det == "LSD"
        Pn = Sn_LSD(f) / 5.0
    elseif det =="NSEMO"
        Pn = Sn_nsemo(f) / 5.0
    else
        detector0 = LIGO(name=det, TObs=TObs);
        Pn = detector0.Sn(f) / 5.0
    end
    Ωhs = interpolate((fres(ms),), ΩgwPBHlog.(1e-3 , fres(ms), ms , mc, σ, "old"), Gridded(Linear()))
    if 1.0*fres(ms)[1] < 1.0*f < 1.0*fres(ms)[end]
        if det == "LIGO"
            Rf = γLIGO(f)/5  #R_eff(f) = Gamma(f)
        else
            Rf = 1.0/5   #Gaamma(f) = 1
        end
        Ωh = Ωhs(f)
        Sh = 3.0*H0_SI^2.0 * Ωh / 2.0 /f^3.0 / π^2.0
        num = Sh^2.0 * Rf^2.0
        den = (1.0/25+Rf^2.0)*Sh^2.0+Pn^2.0 + 2.0*Sh*Pn/5.0
        Tobs = TObs*365.0*86400.0
        return num/den*Tobs
    else
        return 0
    end
end

@memoize function SNRPBHlog_ms(ms  ,  mc, σ , det , TObs = 4.0)
    if det == "LSD"
        fmin = 1.0e4
        fmax = 3.0e5
    elseif det =="NSEMO"
        fmin = h_nx[1]
        fmax = h_nx[end]
    else
        detector = LIGO(name=det, TObs=TObs);
        fmin = detector.fPlotRange[1]
        fmax = detector.fPlotRange[2]
    end
    
    SNRPBHlog_ms_tobeint(f, ms , mc, σ  , det , TObs) = (fmax-fmin) * SNRPBHlog_int0(fmin + (fmax-fmin)*f , ms  , mc, σ  ,det , TObs)

    function integrand(x, f)
        f[1] = SNRPBHlog_ms_tobeint(x[1], ms , mc, σ  , det , TObs)
    end
    result, err = cuhre(integrand ,  rtol=1e-2 , maxevals=Int(1e9))
    result[1]^(0.5)
end
######################################################################################################################################

function smoothS(M)
    mmin = 3.96
    δm   = 4.83    
    if M<mmin
        return 0
    elseif M >= (mmin + δm)
        return 1
    else
        m0 = M-mmin
        fterm = exp(δm/m0+δm/(m0-δm))
        return (fterm+1)^(-1)
    end
end

function Rbirth_int0(z,M,M2,Td)
    if M>=M2
        td=exp(1.0*Td)
        t = ztot_inter0(z)
        z1 = ttoz_inter0(t-td)
        num = ndot_ligo_inter(M,z1)
        mmin = 3.96
        δm   = 4.83  
        mmax=87.14
        b=0.43
        mbreak = mmin +b*(mmax-mmin)
        
        if mmin<M<mbreak
            pM = M^(-1.58)*smoothS(M) #α1=1.58
        elseif mbreak<M<mmax
            pM = M^(-5.59)*smoothS(M) #α2=5.59
        else
            pM = 0
        end
        
        if pM==0
            return 0
        else
            βq = 1.40
            pm2 = (M2/M)^βq*smoothS(M2)
            return num*pm2*pM
        end
    else
        return 0
    end
end

function Rbirth_unnorm(z,M,M2)
        Tmin = log(50.0e6)
        Tmax = log(1/H0)
        Rbirth_tobeint(z,M,M2, Td) = Rbirth_int0(z,M,M2, Tmin + (Tmax-Tmin)*Td)* (Tmax-Tmin)
        function integrand(x, f)
            f[1] = Rbirth_tobeint(z,M, M2, x[1])
        end
        result, err = cuhre(integrand ,  rtol=1e-1 , maxevals=Int(1e9))
        result[1]
end

function Rlocal_unnorm1(z,M)
    Mmin=3.0
    Mmax=M
#     Rlocal_tobeint(z,M,M2) = (Mmax-Mmin)^2 * Rbirth_unnorm(z,Mmin + (Mmax-Mmin)*M, Mmin + (Mmax-Mmin)*M2)
    Rlocal_tobeint(z,M,M2) = (Mmax-Mmin) * Rbirth_unnorm(z,M, Mmin + (Mmax-Mmin)*M2)
    function integrand(x, f)
        f[1] = Rlocal_tobeint(z,M,x[1])
    end
    result, err = cuhre(integrand ,  rtol=1e-1 , maxevals=Int(1e9))
    result[1]
end

@memoize function Rlocal_unnorm(z)
    Mmin=3.01
    Mmax=50.0
#     Rlocal_tobeint(z,M,M2) = (Mmax-Mmin)^2 * Rbirth_unnorm(z,Mmin + (Mmax-Mmin)*M, Mmin + (Mmax-Mmin)*M2)
    Rlocal_tobeint2(z,M) = (Mmax-Mmin) * Rlocal_unnorm1(z,Mmin + (Mmax-Mmin)*M)
    function integrand(x, f)
        f[1] = Rlocal_tobeint2(z,x[1])
    end
    result, err = cuhre(integrand ,  rtol=1e-1 , maxevals=Int(1e9))
    result[1]
end

###  Rlocal_unnorm(0.0)=8.218350200746155e-6 for P(m1)\propto m^(-2.35)
###  Rlocal_unnorm(0.0)=0.004039008480475333 for p(m1),p(m2) uniform distribution
###  Rlocal_unnorm(0.0)=1.5966933724289513e-5 for broken pl
Rlocal_sBH=1.5966933724289513e-5

@memoize Rlocal(z)=Rlocal_unnorm(z)/Rlocal_sBH*23.9e-9 #yr^-1 Mpc^-3
# @time Rlocal(2.0)

@memoize Rbirth(z,M1,M2) = Rbirth_unnorm(z,M1,M2)/Rlocal_sBH*23.9e-9
# Rbirth(z,M)=Rbirth_unnorm(z,M)/Rbirth_unnorm(0,M)*56
