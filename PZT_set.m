function PZT_Set(z0, z1, dt);

%  PI-PiezoTranslationController: sets z_position to z1 (in microns)

% (c) Ingo Gregor 2004

% Defininitions


ao = analogoutput('nidaq',2);            %  Get HW access
chans = addchannel(ao,0);                %  open an output channel

set(ao,'SampleRate',1000)                %  set sample rate for output values 
rate = get(ao,'SampleRate');             %  check it 

                                         %  create series of output values
                                         
data = 1 + sin(linspace(-pi/2, pi/2, rate*dt))';    
data = (z1-z0)/20 * data + z0/10;

putdata(ao,data);                        %  store in outpuf buffer
start(ao);                               %  send trigger to start D/A conversion

waittilstop(ao,3*dt);                    %  wait until finished 

delete(ao);                              % clean up
clear all;