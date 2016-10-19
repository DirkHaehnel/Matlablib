% -------Warren D. Smith 1998-------------
% This useful file causes MAPLE to know about
% fundamental physical constants, useful data, and units.

% I feel cgs, emu, esu units are an abomination.
% If you want to have no units, go totally dimensionless
% (for example, the Planck units). If you want to have units,
% use SI. The cgs, emu, esu systems are a silly halfassed attempt
% to have it both ways and get rid of some units and keep others.

%    The present file features conversion factors to random non-SI
% units like foot, horsepower, pound, ounce. The entire group
% of units is known to it, so all quantities are dimensioned
% automatically. Meter, second, coulomb, kelvin, and kg are
% my fundamental units. (Arguably kelvin is not needed.
% Also for radiation dose units, I use the "disintegration", e.g. 
% disintegration/second.)

%   I have so far been unable to stomach the "candela".
% The old definition was...
% candle =  1/60 the luminous intensity of 1 cm^2 of a Blackbody at the
% melting Temperature of platinum (2042.15 kelvin) = 1.6437 *watt
% = meter^2 / 600000 * sigmasb * platinummelttemp^4;
% a candela is the light intensity from a candle (now redefined
% as "1/683 watt per steradian" of 555nm wavelength light (which the eye
% is maximally sensitive to). 
% Note, this seems to be a lot smaller than the old defn.
% Huh? Then   lumen = candela*steradian = 1/683 *watt;  
% lux = lumen/meter^2 = 1/683 *watt/meter^2; candlepower = candela.
% Yuk!

% I also hate the "radian" and "steradian"; they should be
% dimensionless, in my view.

%    Finally, decibels and bels are another abomination.
%   decibels = 10*bels = 10*log10(Intensity/ReferenceIntensity);
% and for use in auditory work, one uses
% ReferenceIntensity = 1 phon = 10^(-16) *watt/cm^2/sec at 1000 *hertz.
% Idiotic - haven't these guys heard of scientific notation to deal with
% logarithms?

%   This MAPLE code does not know about precision (usually the
% last 2 digits of our numbers are imprecise). My values
% are CODATA values plus miscellaneous sources like
% CRC handbook. If you really care about errors you have to know
% about the correlations among the errors in the 
% fundamental constants (many errors are highly correlated)
% and that would require a lot of extra software, which I
% didn't feel like trying to write.

% Update 1999; 
% PJ Mohr BN Taylor
% CODATA recommended values of the fundamental physical constants 1998
% J Phys Chem Ref Data 28 (1999) 
% finally came out
% http;//physics.nist.gov/cgi-bin/cuu/
%   http;//physics.nist.gov/PhysRefData/contents.html

% The  1998  CODATA
% Recommended Values of the Fundamental Physical
% Constants, Web Version 3.0," available at
% http;//physics.nist.gov/constants (National Institute of
% Standards and Technology, Gaithersburg, MD 20899,
% release date 23 July 1999).
% http;//physics.nist.gov/cuu/Constants/index.html
% http;//physics.nist.gov/constants
% http;//www.codata.org/ 
% http;//pdg.lbl.gov/1998/astrorpp.ps
% these sites very recommended - can now download any correlation coeff too

% A detailed description of the data and its analysis that led to these
% values will be published by early 2000 in the Journal of Physical
% and Chemical Reference Data.

% "Review of particle properties", Phys Rev D50 (1994) 1-.
% PRD 54 (July 1996) 1-720.
% C.Caso et al.; Review of particle physics, Europ. Phys. J. C3 (1998) 1-.
% http;//pdg.lbl.gov

% ---------units (and data);---------------

meter = 1;
second = 1;
coulomb = 1; 
kelvin = 1;
kg = 1;
disintegration = 1;
radian = 1;
steradian = 1;

% prefixes;
zepto = 10^(-21);
atto = 10^(-18);
femto = 10^(-15);
pico = 10^(-12);
nano = 10^(-9);
micro = 10^(-6);
milli = 10^(-3);
centi = 10^(-2);
kilo = 10^3;
mega = 10^6;
giga = 10^9;
tera = 10^12;
peta = 10^15;
exa = 10^18;
zetta = 10^21;

joule = meter^2 * kg / second^2;
watt = joule/second;
mole = 6.02214199 * 10^(23);  % 79 ppb uncert
ampere = coulomb/second;
hertz = 1/second;
newton = meter*kg/second^2;   % force
pascal = newton/meter^2;      % pressure
poise = pascal*second;        % absolute viscosity
stokes = meter^2 /second;     % kinematic viscosity
volt = watt/ampere;
farad = coulomb/volt;
ohm = volt/ampere;
mho = 1/ohm;
siemens = 1/ohm;
weber = volt*second;
tesla = weber/meter^2;
henry = weber/ampere;

revolution = 2*pi*radian;
angulardegree = pi/180*radian;
angularminute = angulardegree/60;
angularsecond = angularminute/60;

sipole = meter*coulomb/second; % my own unit of magentic monopole charge
% One SIpole is the amount of magnetic charge that would
% experience a force of 1 Newton if placed in a 1 Tesla magnetic
% field.

becquerel = disintegration/second;  % radioactivity unit
gray = joule/kg;  % for measuring absorbed radiation dose
sievert = gray;
% a sievert is the amount of radiation equivalent in health effect to 1 gray
% of some kind of standard type of radiation.

minute = 60*second;
hour = 60*minute;
solarday = 24*hour;  % increases 1 millisecond/century
starday = 86164.09 * second;
year = 365.24219879  *solarday;
month = year/12;
century = 100*year;

denier = kg/meter;
dyne = 10^(-5) *newton;
erg = 10^(-7) *joule;

gauss = 10^(-4) * tesla;
maxwell = 10^(-8) *weber;

horsepower = 746 * watt;  % several definitions exist 735 to 746 watt.

jansky = 10^(-26) * watt * second / meter^2;

oersted = 7.957747 * ampere/meter; % says CRC handbook

lb = .45359237  *kg;  % exact
poundmass = lb;
ounce = lb/16;
shortton = 2000 *lb;
longton =  2240 *lb;
grain = lb/7000;
slug = 14.59390 *kg;

lbf = 4.4482216152605 * newton; % exact
poundforce = lbf;
torr = 133.322 *pascal;
mmhg = torr;

liter = 10^(-3) *meter^3;

sverdrup = 10^6 *meter^3/second; % volume flux. Used in oceanography.
kayser = 1/(centi*meter); % used in spectroscopy
cc = (centi*meter)^3;
micron = 10^(-6) * meter;

metricton = 10^3 *kg;  % also "tonne"
% electron-volt
ev = 1.602176462 * 10^(-19) *joule;
km = 10^3 *meter;
mm = 10^(-3) *meter;
cm = 10^(-2) *meter;
astronomicalunit = 149597900 *km;  % earth-sun distance
parsec = cot(1 * angularsecond) * astronomicalunit; % 3.085678*10^(16) *meter; 3.26*lightyear
nauticalmile = 1852*meter;
% nauticalmile = 6076.115486*foot;

angstrom = 10^(-10) * meter;
barn = 10^(-28) *meter^2;
fermi = 10^(-15) *meter;

bar = 10^(5) *pascal;

abampere = ampere/10;
abcoulomb = coulomb/10;
statampere = 3.335641 * 10^(-10) *ampere;
statcoulomb = 3.335641 * 10^(-10) *coulomb;
curie = 3.7*10^(10) *becquerel;
roentgen = 2.58 * 10^(-4) *coulomb/kg;
rad = 10^(-2) *gray;
remunit = 10^(-2) *sievert;

fullsphere = 4*pi*steradian;

knot = nauticalmile/hour;
hectare = 10^4 *meter^2;
angstrom = 10^(-10)*meter;
inch = 2.54 * 10^(-2)*meter;
foot = 12*inch;
yard = 3*foot;
printerpoint = 0.13837 *inch;  % or is it inch/72? close
mile = 5280*foot;
acre = mile^(2) /640; % acre = 43560 *foot^2

psi = lbf/inch^2;  % 6895*pascal

usgallon = 231 * inch^3;
brgallon = 277.420 *inch^3;
quart = usgallon/4;
pint = quart/2;
cup = pint/2;
gill = cup/2; % british cups, pints, gills are each a bit bigger
fluidounce = cup/8;
tablespoon = fluidounce/2;
teaspoon  =  tablespoon/3;
boardfoot = 144 *inch^3;

gram = 10^(-3) *kg;
carat = .2 *gram;

btu = 1.055056 * 10^(3) * joule;
ustherm = 1.054804 * 10^(8) *joule;  % used by US natural gas industry

watertripletemp = 273.16 *kelvin; % exact

watertriplepressure = 611.73 *pascal;
platinummelttemp = (1769.0+273.15) *kelvin;
waterfreezetemp = 273.15 *kelvin;
waterboiltemp = 373.15 *kelvin;
% celsius = waterfreeze + kelvin. fahrenheit = celsius * (5/9) + 32;

caloriethermo = 4.184 *joule;
calorieintl = 4.1868 *joule;
calorie = calorieintl;

% weapons kiloton is by convention 10^(12) calories.
kiloton = 10^(12) * calorie;
mageaton = 10^(3) * kiloton;

rankine = (5/9) *kelvin;
standardatmosphere = 101325 *pascal;
standardgravity = 9.80665 *meter/second^2;
vsoundair = 331.4 *meter/second;
specheatair = 1010 *joule/kg/kelvin; % at const pressure
densityair = 1.21 * kg/meter^3;
vsoundwater = 1496.7 *meter/second;
densitywater = 1000*kg/meter^3;
specheatwater = 4190 *joule/kg/kelvin; % at const pressure
heatfusionwater = 333 * kilo*joule/kg;
heatvapwater = 2260* kilo*joule/kg;
indexrefwater = 1.33; % at 589 nm

resistivitysilver = 1.587 * 10^(-8) * ohm * meter; % at 293 Kelvin

% ------------------fundamental constants;-----------------

% these constants suffice to define the planck units;

% speed of light in vacuum    exact
c = 299792458 * meter/second;
lightyear = c*year;

% boltzmann constant
kb = 1.3806503 * 10^(-23) *joule/kelvin;  % 1700 ppb uncert

% permittivity &amp; permeability of vacuum
mu0 = 4*pi*1.0 * 10^(-7) *henry/meter; % exact
eps0 = 1/(c^2 * mu0); % exact
% eps0 = 8.854187817 * 10^(-12) *farad/meter
% mu0*eps0*c^2=1

% newton gravitational constant (force = G*m1*m2/r^2)
G = 6.673 * 10^(-11) *newton*meter^2/kg^2;  
% +- .010; huge relative uncert .0015

% coulomb law constant (force = coulombconst * q1*q2/r^2)
coulombconst = 1/ (4 * pi *eps0);

% planck's constant
h = 6.62606876 * 10^(-34) *joule*second;
hbar = 1.054571596 * 10^(-34) *joule*second;
% hbar = h/(2*pi);
% 78 ppb uncert

% ------------planck length, mass, time, charge;--------------
% These are the fundamental units of 
% length, mass, time, charge, temperature
% in a unit system where c=G=hbar=kb=1 and 
% where coulombconst*qplanck^2 = G*mplanck^2.

lplanck = sqrt(hbar*G/c^3);
mplanck = sqrt(hbar*c/G);
tplanck = sqrt(hbar*G/c^5);
qplanck = sqrt(4*pi*hbar*c*eps0);
% qplanck / qe = 1/sqrt(alpha), reciprocal sqrt of fine structure constant!

tempplanck = (mplanck*c^2)/kb;

%  lplanck    = 1.616048609*10^(-35) *meter;
%  tplanck    = 5.390557920*10^(-44) *second;
%  mplanck    = 2.176714074*10^(-8)  *kg;
%  qplanck    = 1.875546788*10^(-18) *coulomb; 
%  tempplanck = 1.416957021*10^(32)  *kelvin;

%  qplanck = 11.7062 charge quanta = sqrt(137).
%  mplanck = 2.389527889 * 10^(22) electron masses.

% ------other fundamental constants;-----------------

% charge of proton
qe = 1.602176462 * 10^(-19) *coulomb;
% 39 ppb uncert

% atomic mass unit
amu = 10^(-3) * kg/mole;
% amu = 1.66053873*10^(-27) *kg;  % 79 ppb uncert
% mole is number of carbon-12 atoms in 12 grams, amu is
% 1/12 mass of a carbon-12 atom.

% rest mass of electron, proton, neutron, muon, deuteron, alpha particle
% me = 9.10938188 * 10^(-31) *kg;   % 79 ppb uncert
me = 5.48579903 * 10^(-4) * amu;   %  +-   .00000013
me = .510998902 * mega * ev/c^2;  % 40 ppb uncert
% failed excited lepton searches; e &gt;86 GeV, muon &gt;=86 GeV, tau&gt;72 GeV

mp = 1.67262158* 10^(-27) *kg;  % 79 ppb uncert
mp = 938.271998 * mega * ev / c^2;   % 40 ppb uncert
mp = 1.00727646688 * amu;  %  .13ppb uncert

mpbyme = 1836.1526675; % 2.1ppb uncert

% mn = 1.6749286 * 10^(-27) *kg;
mn = 1.008664904 * amu;
%     +-.000000014

% mmu = 1.8835327*10^(-28) *kg;
mmu = .113428913 * amu;  %  +-   .000000017
mmu = 105.658389 * mega * ev/c^2; % +-.000034  (1998 value)

mdeuteron = 3.3435860 * 10^(-27) *kg;
mdeuteron = 1875.612762 * mega * ev/c^2; % 40 ppb
mtriton   = 3.015500 * amu;
m3he      = 3.014932 * amu;
malpha = me / (1.37093354 * 10^(-4));
% malpha =       4.001505 * amu;
mpi  = 135.0 * mega*ev/c^2; % mass of pi
mpi0 = 139.6 * mega*ev/c^2; % mass of pi_0^{+-}
mtau = 1777.05 * mega * ev / c^2; % +.29, -.26  (1998 value)

% heavy lepton searches supposedly go up to 81.5GeV now with 95% conf 
% nothing charged there. Ackerstaff Euro Phys J C1 (1998) 45 
% also 95% conf nothing uncharged below 69GeV

mquarku = 3.25 * mega * ev/c^2; % charge +2/3  (1.5 to 5)
mquarkd = 6    * mega * ev/c^2; % charge -1/3 % ??? wrong way ???  (3 to 9)
mquarks = 115 * mega * ev/c^2; % charge -1/3  (60 to 170)
mquarkc = 1.25  * giga * ev/c^2; % charge +2/3  (1.25 +-.15)
mquarkb = 4.25  * giga * ev/c^2; % charge -1/3  (4.25 +- .15)
mquarkt = 174.3 * giga * ev/c^2; % charge +2/3. 174.3+-5.1
% all quarks spin 1/2, baryon %  1/3, RGB or anti-RGB color.

mw = 80.41  * giga * ev/c^2; % +-.4 Gev. Charge +-1. Spin 1. W boson.
mz = 91.187 * giga * ev/c^2; % +-.031 Gev. Charge 0. Spin 1. Z boson.

% particle mean lifetimes
meanlifemu = 2.19703 * micro * second;% uncert .00004
meanlifetau = .2900 * pico * second;  % +-.0012
meanlifen = 888.6 * second; % +- 3.5 second
% meanlifep &gt; 10^31 years for some investigated modes

% Bohr magneton
bohrmag = qe*hbar/(2*me);
% bohrmag = mue/1.001159652;
% bohrmag = 9.2740154 * 10^(-24) *joule/tesla;

% Nuclear magneton
nucmag = qe*hbar/(2*mp);;
% nucmag = 5.0507866*10^(-27) * joule/tesla;

% magnetic moment of electron, proton, neutron, muon, deuteron
% mue = 9.2847701 * 10^(-24) *joule/tesla;
mue = 1.001159652193 * bohrmag;
%  +-    .000000000010

% mup = 1.41060761 * 10^(-26) *joule/tesla;
mup = 2.792847386 * nucmag;
%  +-    .000000063

% mun = 9.6623707 * 10^(-27) *joule/tesla;
mun = 1.91304275 * nucmag; % has negative sign, incidentally
%      +-.00000045

mumu = 4.4904514 * 10^(-26) *joule/tesla;
mudeuteron = 4.3307375 * 10^(-27) *joule/tesla;

% ------------"atomic units";--------------
% mass     me
% charge   qe
% action   hbar
% length   rbohr
% energy   hartree

% Bohr radius (estimate, for hydrogen atom)
rbohr = 4*pi*eps0*hbar^2/(me*qe^2);
% rbohr = 5.29177208 * 10^(-11) *meter; % 3.7ppb

% hartree energy unit, about 27 ev
% hartree = 2*rydberg*h*c = 4.3598 * 10^(-18) *joule;
hartree = hbar^2 / (me*rbohr^2);

atomictimeunit = hbar / hartree;
atomicvelocityunit = rbohr / atomictimeunit;   % about c/137
atomicforceunit = hartree / rbohr;
atomicmomentumunit = hbar / rbohr;
atomiccurrentunit = qe / atomictimeunit;
atomicefieldunit = hartree/(qe*rbohr);
atomicpotentialunit = hartree/qe;
atomicedipoleunit = qe*rbohr;
atomicmagfieldunit = hbar / (qe*rbohr^2);

atomicmagdipoleunit =  2*bohrmag;

% ------derived constants;-----------------

% magnetic flux quantum
phi0 = h/(2*qe);
% phi0 = 2.06783461 * 10^(-15) *weber

% alpha = fine structure constant
% alpha = 1/137.0359895;
alpha = qe^2 / (4*pi*eps0*hbar*c);

% stefan-boltzmann constant
sigmasb = pi^2*kb^4/(60*hbar^3*c^2);
% sigmasb =  5.670400*10^(-8) *watt/meter^2/kelvin^4; % 7000 ppb

% first radiation constant
c1radiation = 2*pi*h*c^2;
% c1radiation = 3.7417749 * 10^(-16) watt/meter^2;
c2radiation = h*c/kb;
% c2radiation = 0.01438769 *meter*kelvin;
% The max intensity lambda in blackbody radiation is
% wienconst*T*cm/kelvin, where
wienconst = c2radiation / 4.965114231;
% wienconst = 0.2897756 *cm/kelvin;

% rydberg (approx hydrogen waves per meter, infinite mass nuc)
rydberg = me*c*alpha^2/(2*h);
% rydberg = 10973731.56834(24) /meter;

% classical electron radius
% re = 2.817940285 * 10^(-15) *meter; % 11 ppb
re = alpha^2 * rbohr;

% compton wavelengths of electron, proton, neutron, muon
e_compton = h/(me*c);
p_compton = h/(mp*c);
n_compton = h/(mn*c);
mu_compton = h/(mmu*c);

fermicouplingconst = 1.16639 * 10^(-5) / (giga * ev)^2; % 9000 ppb
sinsquaredweakmixingangle = 0.23124;  % .001 relative uncert
% codata says .2235 +- .0023

strongcouplingconst = .119;

% Astronomy

radgeosync = 42200 * kilo * meter;
vescapeearth = 11200 * meter/second;
muearth = 8.0*10^(22) *joule/tesla; % magnetic dipole moment
efieldearth = 150 * volt/meter;  % down. At earth's surface. mean.

suntemp = 5780*kelvin;
sunluminosity = 3.846 * 10^(26) * joule/second;
solarconst = sunluminosity / (4*pi*astronomicalunit^2); % varies 2%

mmilkyway = 2.2*10^(41) *kg;
moceans = 1.4 *10^(21) *kg;

velocorbitalgalaxy = 220 * km/second; % of solar system round milky way galaxy
velocwrtcosmicbackground = 369.3 * km/second; % +-2.5

disttocenterofmilkyway = 8.0 * mega * parsec;

distproximacentauri = 4.22 * 299792458 * year;
% **Alpha  Centauri system;
% Separation between  Alpha  Centauri A and B varies from 11 to 35 AU;
% they take 80 years to orbit around each other.
% Proxima  Centauri is currently 13,000 AUs from A and B; no orbital
% parameters of it are known.
%   A            B           proxima     sun
%  5800K        5300K       2700K      5800K
%   4.35         4.35        4.22         0 lightyears away
%   Yellow       Orange         Red    Yellow
distandromedagalaxy = 2*10^(22) *meter;
distlargemagellaniccloud = 55*kilo*parsec;
distvirgocluster = 20*mega*parsec;
disttouniverseedge = 10^(26) * meter;
distearthmoon = 384404 *km;  % +3 cm per year
periodmoon = 27 *solarday;  % ??? +2 milliseconds/century 

distsunmercury = .38  * astronomicalunit;
distsunvenus   = .72  * astronomicalunit;
distsunearth   = 1.0  * astronomicalunit;
distsunmars    = 1.52 * astronomicalunit;
distsunjupiter = 5.20 * astronomicalunit;
distsunsaturn  = 9.5  * astronomicalunit;
distsunuranus  = 19.2 * astronomicalunit;
distsunneptune = 30.0 * astronomicalunit;
distsunpluto   = 39.5 * astronomicalunit;
% pluto has 1 satellite "charon" with about .1 its mass; 
% it has day 6.38 earth days.

% radiusearth   = 6371.006 *km;
radiusearth   = 6378.140 *km;% equatorial
radiusmoon    = 1738 *km;

radiusmercury = .38  * radiusearth;
radiusvenus   = .96  * radiusearth;
radiusmars    = .53  * radiusearth;
radiusjupiter = 10.8 * radiusearth;
radiussaturn  = 9.0  * radiusearth;
radiusuranus  = 4.1  * radiusearth;
radiusneptune = 3.85 * radiusearth;
radiuspluto   = 1137 *km;  % +-8
radiussun     = 695990 *km;  % equatorial
radiussunschwarzschild     = 2.95325008 * km;
radiusmilkyway = 15 * kilo *parsec;
thicknessmilkyway = 1 * kilo *parsec;

mmercury = .33022 * 10^(24) * kg;
mvenus   = 4.8690 * 10^(24) * kg;
mearth   = 5.97370 * 10^(24) * kg;
mmoon    = .073483 * 10^(24) * kg;
mmars    = .64191 * 10^(24) * kg;
mjupiter = 1898.8 * 10^(24) * kg;
msaturn  = 568.50 * 10^(24) * kg;
muranus  = 86.625 * 10^(24) * kg;
mneptune = 102.78 * 10^(24) * kg;
mpluto   = 1.27 * 10^(22) * kg;
msun = 1.98892 * 10^(30) *kg;

hubnum = .5; % could be anywhere .2 to 1.0
hubbleconst = 70*hubnum*km/second/(mega*parsec);
% http;//ucsu.colorado.edu/~lisle/main.html
% ``The cosmological implications of Hipparcus''
% direct trig parallax distance measurement to nearby Cepheids with
% Hipparcus satellite. Leads to 10% revision of Hubble constant, now
% about 60  km/s/Mpc.
% another place claims 60-80 is the range; I'm going with it.

criticaldensity = 3*hubbleconst^2/(8*pi*G); % needed to close the universe
% matter density is supposed to be .2 to 1.0 times this

entropydensityofuniverse = 2899.3 * (2.853/2.728)^3 * kb / cc; 

tempuniversebackground = 2.853 * kelvin;  % plus or minus .002

% various rough estimates of universe-quantities on Misner-Thorne-Wheeler p738;
% max radius 19*10^9 lightyear
% time start to max 3*10^10 year
% time bang to crunch 6*10^10 year
% age now 10^10 year
% radius now 13*10^9 lightyear    note; larger than age
% volume now 3.83*10^79 meter^3
% density at maximum   5*10^(-30) * gram/cc
% rate of increase of radius now   .66 lightyear/year
% amount of matter   5.68*10^(53) kg
% baryon number      3.39*10^(80)

% MTW 29.6p796
% density of universe now is &lt;10^(-28) g/cc and &gt;2*10^(-31) g/cc

% MTW when universe about 10^5 years old, density was 10^(-20) g/cc and temp 3000K
% and then; 
% 1. radiation ceased to domaine, now matter
% 2. universe became transparent
% 3. hydrogen atoms began to form (no longer ionized perpetually)

