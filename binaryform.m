% artificial binary operation 
function y=binaryform(x)
y=zeros(1,8);
if(x>127), y(8)=1; x=x-128; end;    
if(x>63), y(7)=1; x=x-64; end;    
if(x>31), y(6)=1; x=x-32; end;    
if(x>15), y(5)=1; x=x-16; end;    
if(x>7), y(4)=1; x=x-8; end;    
if(x>3), y(3)=1; x=x-4; end;    
if(x>1), y(2)=1; x=x-2; end;    
if(x>0), y(1)=1;  end;    