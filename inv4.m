function R = inv4(M,N)

R(:,1,1) = (-M(:,2,2).*M(:,3,4).*M(:,4,3).*N(:,1,1)+M(:,2,2).*M(:,3,3).*M(:,4,4).*N(:,1,1)+M(:,1,4).*M(:,3,3).*M(:,4,2).*N(:,2,1)-M(:,1,3).*M(:,3,4).*M(:,4,2).*N(:,2,1)-M(:,1,4).*M(:,3,2).*M(:,4,3).*N(:,2,1)+M(:,1,2).*M(:,3,4).*M(:,4,3).*N(:,2,1)+M(:,1,3).*M(:,3,2).*M(:,4,4).*N(:,2,1)-M(:,1,2).*M(:,3,3).*M(:,4,4).*N(:,2,1)+M(:,1,4).*M(:,2,2).*M(:,4,3).*N(:,3,1)-M(:,1,3).*M(:,2,2).*M(:,4,4).*N(:,3,1)-M(:,1,4).*M(:,2,2).*M(:,3,3).*N(:,4,1)+M(:,1,3).*M(:,2,2).*M(:,3,4).*N(:,4,1)+M(:,2,4).*((M(:,1,3).*M(:,4,2)-M(:,1,2).*M(:,4,3)).*N(:,3,1)+M(:,3,3).*(-M(:,4,2).*N(:,1,1)+M(:,1,2).*N(:,4,1))+M(:,3,2).*(M(:,4,3).*N(:,1,1)-M(:,1,3).*N(:,4,1)))+M(:,2,3).*((-M(:,1,4).*M(:,4,2)+M(:,1,2).*M(:,4,4)).*N(:,3,1)+M(:,3,4).*(M(:,4,2).*N(:,1,1)-M(:,1,2).*N(:,4,1))+M(:,3,2).*(-M(:,4,4).*N(:,1,1)+M(:,1,4).*N(:,4,1))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,1,2) =(-M(:,2,2).*M(:,3,4).*M(:,4,3).*N(:,1,2)+M(:,2,2).*M(:,3,3).*M(:,4,4).*N(:,1,2)+M(:,1,4).*M(:,3,3).*M(:,4,2).*N(:,2,2)-M(:,1,3).*M(:,3,4).*M(:,4,2).*N(:,2,2)-M(:,1,4).*M(:,3,2).*M(:,4,3).*N(:,2,2)+M(:,1,2).*M(:,3,4).*M(:,4,3).*N(:,2,2)+M(:,1,3).*M(:,3,2).*M(:,4,4).*N(:,2,2)-M(:,1,2).*M(:,3,3).*M(:,4,4).*N(:,2,2)+M(:,1,4).*M(:,2,2).*M(:,4,3).*N(:,3,2)-M(:,1,3).*M(:,2,2).*M(:,4,4).*N(:,3,2)-M(:,1,4).*M(:,2,2).*M(:,3,3).*N(:,4,2)+M(:,1,3).*M(:,2,2).*M(:,3,4).*N(:,4,2)+M(:,2,4).*((M(:,1,3).*M(:,4,2)-M(:,1,2).*M(:,4,3)).*N(:,3,2)+M(:,3,3).*(-M(:,4,2).*N(:,1,2)+M(:,1,2).*N(:,4,2))+M(:,3,2).*(M(:,4,3).*N(:,1,2)-M(:,1,3).*N(:,4,2)))+M(:,2,3).*((-M(:,1,4).*M(:,4,2)+M(:,1,2).*M(:,4,4)).*N(:,3,2)+M(:,3,4).*(M(:,4,2).*N(:,1,2)-M(:,1,2).*N(:,4,2))+M(:,3,2).*(-M(:,4,4).*N(:,1,2)+M(:,1,4).*N(:,4,2))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,1,3) =(-M(:,2,2).*M(:,3,4).*M(:,4,3).*N(:,1,3)+M(:,2,2).*M(:,3,3).*M(:,4,4).*N(:,1,3)+M(:,1,4).*M(:,3,3).*M(:,4,2).*N(:,2,3)-M(:,1,3).*M(:,3,4).*M(:,4,2).*N(:,2,3)-M(:,1,4).*M(:,3,2).*M(:,4,3).*N(:,2,3)+M(:,1,2).*M(:,3,4).*M(:,4,3).*N(:,2,3)+M(:,1,3).*M(:,3,2).*M(:,4,4).*N(:,2,3)-M(:,1,2).*M(:,3,3).*M(:,4,4).*N(:,2,3)+M(:,1,4).*M(:,2,2).*M(:,4,3).*N(:,3,3)-M(:,1,3).*M(:,2,2).*M(:,4,4).*N(:,3,3)-M(:,1,4).*M(:,2,2).*M(:,3,3).*N(:,4,3)+M(:,1,3).*M(:,2,2).*M(:,3,4).*N(:,4,3)+M(:,2,4).*((M(:,1,3).*M(:,4,2)-M(:,1,2).*M(:,4,3)).*N(:,3,3)+M(:,3,3).*(-M(:,4,2).*N(:,1,3)+M(:,1,2).*N(:,4,3))+M(:,3,2).*(M(:,4,3).*N(:,1,3)-M(:,1,3).*N(:,4,3)))+M(:,2,3).*((-M(:,1,4).*M(:,4,2)+M(:,1,2).*M(:,4,4)).*N(:,3,3)+M(:,3,4).*(M(:,4,2).*N(:,1,3)-M(:,1,2).*N(:,4,3))+M(:,3,2).*(-M(:,4,4).*N(:,1,3)+M(:,1,4).*N(:,4,3))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,1,4) =(-M(:,2,2).*M(:,3,4).*M(:,4,3).*N(:,1,4)+M(:,2,2).*M(:,3,3).*M(:,4,4).*N(:,1,4)+M(:,1,4).*M(:,3,3).*M(:,4,2).*N(:,2,4)-M(:,1,3).*M(:,3,4).*M(:,4,2).*N(:,2,4)-M(:,1,4).*M(:,3,2).*M(:,4,3).*N(:,2,4)+M(:,1,2).*M(:,3,4).*M(:,4,3).*N(:,2,4)+M(:,1,3).*M(:,3,2).*M(:,4,4).*N(:,2,4)-M(:,1,2).*M(:,3,3).*M(:,4,4).*N(:,2,4)+M(:,1,4).*M(:,2,2).*M(:,4,3).*N(:,3,4)-M(:,1,3).*M(:,2,2).*M(:,4,4).*N(:,3,4)-M(:,1,4).*M(:,2,2).*M(:,3,3).*N(:,4,4)+M(:,1,3).*M(:,2,2).*M(:,3,4).*N(:,4,4)+M(:,2,4).*((M(:,1,3).*M(:,4,2)-M(:,1,2).*M(:,4,3)).*N(:,3,4)+M(:,3,3).*(-M(:,4,2).*N(:,1,4)+M(:,1,2).*N(:,4,4))+M(:,3,2).*(M(:,4,3).*N(:,1,4)-M(:,1,3).*N(:,4,4)))+M(:,2,3).*((-M(:,1,4).*M(:,4,2)+M(:,1,2).*M(:,4,4)).*N(:,3,4)+M(:,3,4).*(M(:,4,2).*N(:,1,4)-M(:,1,2).*N(:,4,4))+M(:,3,2).*(-M(:,4,4).*N(:,1,4)+M(:,1,4).*N(:,4,4))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,2,1) =(M(:,2,1).*M(:,3,4).*M(:,4,3).*N(:,1,1)-M(:,2,1).*M(:,3,3).*M(:,4,4).*N(:,1,1)-M(:,1,4).*M(:,3,3).*M(:,4,1).*N(:,2,1)+M(:,1,3).*M(:,3,4).*M(:,4,1).*N(:,2,1)+M(:,1,4).*M(:,3,1).*M(:,4,3).*N(:,2,1)-M(:,1,1).*M(:,3,4).*M(:,4,3).*N(:,2,1)-M(:,1,3).*M(:,3,1).*M(:,4,4).*N(:,2,1)+M(:,1,1).*M(:,3,3).*M(:,4,4).*N(:,2,1)-M(:,1,4).*M(:,2,1).*M(:,4,3).*N(:,3,1)+M(:,1,3).*M(:,2,1).*M(:,4,4).*N(:,3,1)+M(:,1,4).*M(:,2,1).*M(:,3,3).*N(:,4,1)-M(:,1,3).*M(:,2,1).*M(:,3,4).*N(:,4,1)+M(:,2,4).*((-M(:,1,3).*M(:,4,1)+M(:,1,1).*M(:,4,3)).*N(:,3,1)+M(:,3,3).*(M(:,4,1).*N(:,1,1)-M(:,1,1).*N(:,4,1))+M(:,3,1).*(-M(:,4,3).*N(:,1,1)+M(:,1,3).*N(:,4,1)))+M(:,2,3).*((M(:,1,4).*M(:,4,1)-M(:,1,1).*M(:,4,4)).*N(:,3,1)+M(:,3,4).*(-M(:,4,1).*N(:,1,1)+M(:,1,1).*N(:,4,1))+M(:,3,1).*(M(:,4,4).*N(:,1,1)-M(:,1,4).*N(:,4,1))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,2,2) =(M(:,2,1).*M(:,3,4).*M(:,4,3).*N(:,1,2)-M(:,2,1).*M(:,3,3).*M(:,4,4).*N(:,1,2)-M(:,1,4).*M(:,3,3).*M(:,4,1).*N(:,2,2)+M(:,1,3).*M(:,3,4).*M(:,4,1).*N(:,2,2)+M(:,1,4).*M(:,3,1).*M(:,4,3).*N(:,2,2)-M(:,1,1).*M(:,3,4).*M(:,4,3).*N(:,2,2)-M(:,1,3).*M(:,3,1).*M(:,4,4).*N(:,2,2)+M(:,1,1).*M(:,3,3).*M(:,4,4).*N(:,2,2)-M(:,1,4).*M(:,2,1).*M(:,4,3).*N(:,3,2)+M(:,1,3).*M(:,2,1).*M(:,4,4).*N(:,3,2)+M(:,1,4).*M(:,2,1).*M(:,3,3).*N(:,4,2)-M(:,1,3).*M(:,2,1).*M(:,3,4).*N(:,4,2)+M(:,2,4).*((-M(:,1,3).*M(:,4,1)+M(:,1,1).*M(:,4,3)).*N(:,3,2)+M(:,3,3).*(M(:,4,1).*N(:,1,2)-M(:,1,1).*N(:,4,2))+M(:,3,1).*(-M(:,4,3).*N(:,1,2)+M(:,1,3).*N(:,4,2)))+M(:,2,3).*((M(:,1,4).*M(:,4,1)-M(:,1,1).*M(:,4,4)).*N(:,3,2)+M(:,3,4).*(-M(:,4,1).*N(:,1,2)+M(:,1,1).*N(:,4,2))+M(:,3,1).*(M(:,4,4).*N(:,1,2)-M(:,1,4).*N(:,4,2))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,2,3) =(M(:,2,1).*M(:,3,4).*M(:,4,3).*N(:,1,3)-M(:,2,1).*M(:,3,3).*M(:,4,4).*N(:,1,3)-M(:,1,4).*M(:,3,3).*M(:,4,1).*N(:,2,3)+M(:,1,3).*M(:,3,4).*M(:,4,1).*N(:,2,3)+M(:,1,4).*M(:,3,1).*M(:,4,3).*N(:,2,3)-M(:,1,1).*M(:,3,4).*M(:,4,3).*N(:,2,3)-M(:,1,3).*M(:,3,1).*M(:,4,4).*N(:,2,3)+M(:,1,1).*M(:,3,3).*M(:,4,4).*N(:,2,3)-M(:,1,4).*M(:,2,1).*M(:,4,3).*N(:,3,3)+M(:,1,3).*M(:,2,1).*M(:,4,4).*N(:,3,3)+M(:,1,4).*M(:,2,1).*M(:,3,3).*N(:,4,3)-M(:,1,3).*M(:,2,1).*M(:,3,4).*N(:,4,3)+M(:,2,4).*((-M(:,1,3).*M(:,4,1)+M(:,1,1).*M(:,4,3)).*N(:,3,3)+M(:,3,3).*(M(:,4,1).*N(:,1,3)-M(:,1,1).*N(:,4,3))+M(:,3,1).*(-M(:,4,3).*N(:,1,3)+M(:,1,3).*N(:,4,3)))+M(:,2,3).*((M(:,1,4).*M(:,4,1)-M(:,1,1).*M(:,4,4)).*N(:,3,3)+M(:,3,4).*(-M(:,4,1).*N(:,1,3)+M(:,1,1).*N(:,4,3))+M(:,3,1).*(M(:,4,4).*N(:,1,3)-M(:,1,4).*N(:,4,3))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,2,4) =(M(:,2,1).*M(:,3,4).*M(:,4,3).*N(:,1,4)-M(:,2,1).*M(:,3,3).*M(:,4,4).*N(:,1,4)-M(:,1,4).*M(:,3,3).*M(:,4,1).*N(:,2,4)+M(:,1,3).*M(:,3,4).*M(:,4,1).*N(:,2,4)+M(:,1,4).*M(:,3,1).*M(:,4,3).*N(:,2,4)-M(:,1,1).*M(:,3,4).*M(:,4,3).*N(:,2,4)-M(:,1,3).*M(:,3,1).*M(:,4,4).*N(:,2,4)+M(:,1,1).*M(:,3,3).*M(:,4,4).*N(:,2,4)-M(:,1,4).*M(:,2,1).*M(:,4,3).*N(:,3,4)+M(:,1,3).*M(:,2,1).*M(:,4,4).*N(:,3,4)+M(:,1,4).*M(:,2,1).*M(:,3,3).*N(:,4,4)-M(:,1,3).*M(:,2,1).*M(:,3,4).*N(:,4,4)+M(:,2,4).*((-M(:,1,3).*M(:,4,1)+M(:,1,1).*M(:,4,3)).*N(:,3,4)+M(:,3,3).*(M(:,4,1).*N(:,1,4)-M(:,1,1).*N(:,4,4))+M(:,3,1).*(-M(:,4,3).*N(:,1,4)+M(:,1,3).*N(:,4,4)))+M(:,2,3).*((M(:,1,4).*M(:,4,1)-M(:,1,1).*M(:,4,4)).*N(:,3,4)+M(:,3,4).*(-M(:,4,1).*N(:,1,4)+M(:,1,1).*N(:,4,4))+M(:,3,1).*(M(:,4,4).*N(:,1,4)-M(:,1,4).*N(:,4,4))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,3,1) =(-M(:,2,1).*M(:,3,4).*M(:,4,2).*N(:,1,1)+M(:,2,1).*M(:,3,2).*M(:,4,4).*N(:,1,1)+M(:,1,4).*M(:,3,2).*M(:,4,1).*N(:,2,1)-M(:,1,2).*M(:,3,4).*M(:,4,1).*N(:,2,1)-M(:,1,4).*M(:,3,1).*M(:,4,2).*N(:,2,1)+M(:,1,1).*M(:,3,4).*M(:,4,2).*N(:,2,1)+M(:,1,2).*M(:,3,1).*M(:,4,4).*N(:,2,1)-M(:,1,1).*M(:,3,2).*M(:,4,4).*N(:,2,1)+M(:,1,4).*M(:,2,1).*M(:,4,2).*N(:,3,1)-M(:,1,2).*M(:,2,1).*M(:,4,4).*N(:,3,1)-M(:,1,4).*M(:,2,1).*M(:,3,2).*N(:,4,1)+M(:,1,2).*M(:,2,1).*M(:,3,4).*N(:,4,1)+M(:,2,4).*((M(:,1,2).*M(:,4,1)-M(:,1,1).*M(:,4,2)).*N(:,3,1)+M(:,3,2).*(-M(:,4,1).*N(:,1,1)+M(:,1,1).*N(:,4,1))+M(:,3,1).*(M(:,4,2).*N(:,1,1)-M(:,1,2).*N(:,4,1)))+M(:,2,2).*((-M(:,1,4).*M(:,4,1)+M(:,1,1).*M(:,4,4)).*N(:,3,1)+M(:,3,4).*(M(:,4,1).*N(:,1,1)-M(:,1,1).*N(:,4,1))+M(:,3,1).*(-M(:,4,4).*N(:,1,1)+M(:,1,4).*N(:,4,1))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,3,2) =(-M(:,2,1).*M(:,3,4).*M(:,4,2).*N(:,1,2)+M(:,2,1).*M(:,3,2).*M(:,4,4).*N(:,1,2)+M(:,1,4).*M(:,3,2).*M(:,4,1).*N(:,2,2)-M(:,1,2).*M(:,3,4).*M(:,4,1).*N(:,2,2)-M(:,1,4).*M(:,3,1).*M(:,4,2).*N(:,2,2)+M(:,1,1).*M(:,3,4).*M(:,4,2).*N(:,2,2)+M(:,1,2).*M(:,3,1).*M(:,4,4).*N(:,2,2)-M(:,1,1).*M(:,3,2).*M(:,4,4).*N(:,2,2)+M(:,1,4).*M(:,2,1).*M(:,4,2).*N(:,3,2)-M(:,1,2).*M(:,2,1).*M(:,4,4).*N(:,3,2)-M(:,1,4).*M(:,2,1).*M(:,3,2).*N(:,4,2)+M(:,1,2).*M(:,2,1).*M(:,3,4).*N(:,4,2)+M(:,2,4).*((M(:,1,2).*M(:,4,1)-M(:,1,1).*M(:,4,2)).*N(:,3,2)+M(:,3,2).*(-M(:,4,1).*N(:,1,2)+M(:,1,1).*N(:,4,2))+M(:,3,1).*(M(:,4,2).*N(:,1,2)-M(:,1,2).*N(:,4,2)))+M(:,2,2).*((-M(:,1,4).*M(:,4,1)+M(:,1,1).*M(:,4,4)).*N(:,3,2)+M(:,3,4).*(M(:,4,1).*N(:,1,2)-M(:,1,1).*N(:,4,2))+M(:,3,1).*(-M(:,4,4).*N(:,1,2)+M(:,1,4).*N(:,4,2))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,3,3) =(-M(:,2,1).*M(:,3,4).*M(:,4,2).*N(:,1,3)+M(:,2,1).*M(:,3,2).*M(:,4,4).*N(:,1,3)+M(:,1,4).*M(:,3,2).*M(:,4,1).*N(:,2,3)-M(:,1,2).*M(:,3,4).*M(:,4,1).*N(:,2,3)-M(:,1,4).*M(:,3,1).*M(:,4,2).*N(:,2,3)+M(:,1,1).*M(:,3,4).*M(:,4,2).*N(:,2,3)+M(:,1,2).*M(:,3,1).*M(:,4,4).*N(:,2,3)-M(:,1,1).*M(:,3,2).*M(:,4,4).*N(:,2,3)+M(:,1,4).*M(:,2,1).*M(:,4,2).*N(:,3,3)-M(:,1,2).*M(:,2,1).*M(:,4,4).*N(:,3,3)-M(:,1,4).*M(:,2,1).*M(:,3,2).*N(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*N(:,4,3)+M(:,2,4).*((M(:,1,2).*M(:,4,1)-M(:,1,1).*M(:,4,2)).*N(:,3,3)+M(:,3,2).*(-M(:,4,1).*N(:,1,3)+M(:,1,1).*N(:,4,3))+M(:,3,1).*(M(:,4,2).*N(:,1,3)-M(:,1,2).*N(:,4,3)))+M(:,2,2).*((-M(:,1,4).*M(:,4,1)+M(:,1,1).*M(:,4,4)).*N(:,3,3)+M(:,3,4).*(M(:,4,1).*N(:,1,3)-M(:,1,1).*N(:,4,3))+M(:,3,1).*(-M(:,4,4).*N(:,1,3)+M(:,1,4).*N(:,4,3))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,3,4) =(-M(:,2,1).*M(:,3,4).*M(:,4,2).*N(:,1,4)+M(:,2,1).*M(:,3,2).*M(:,4,4).*N(:,1,4)+M(:,1,4).*M(:,3,2).*M(:,4,1).*N(:,2,4)-M(:,1,2).*M(:,3,4).*M(:,4,1).*N(:,2,4)-M(:,1,4).*M(:,3,1).*M(:,4,2).*N(:,2,4)+M(:,1,1).*M(:,3,4).*M(:,4,2).*N(:,2,4)+M(:,1,2).*M(:,3,1).*M(:,4,4).*N(:,2,4)-M(:,1,1).*M(:,3,2).*M(:,4,4).*N(:,2,4)+M(:,1,4).*M(:,2,1).*M(:,4,2).*N(:,3,4)-M(:,1,2).*M(:,2,1).*M(:,4,4).*N(:,3,4)-M(:,1,4).*M(:,2,1).*M(:,3,2).*N(:,4,4)+M(:,1,2).*M(:,2,1).*M(:,3,4).*N(:,4,4)+M(:,2,4).*((M(:,1,2).*M(:,4,1)-M(:,1,1).*M(:,4,2)).*N(:,3,4)+M(:,3,2).*(-M(:,4,1).*N(:,1,4)+M(:,1,1).*N(:,4,4))+M(:,3,1).*(M(:,4,2).*N(:,1,4)-M(:,1,2).*N(:,4,4)))+M(:,2,2).*((-M(:,1,4).*M(:,4,1)+M(:,1,1).*M(:,4,4)).*N(:,3,4)+M(:,3,4).*(M(:,4,1).*N(:,1,4)-M(:,1,1).*N(:,4,4))+M(:,3,1).*(-M(:,4,4).*N(:,1,4)+M(:,1,4).*N(:,4,4))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,4,1) =(M(:,2,1).*M(:,3,3).*M(:,4,2).*N(:,1,1)-M(:,2,1).*M(:,3,2).*M(:,4,3).*N(:,1,1)-M(:,1,3).*M(:,3,2).*M(:,4,1).*N(:,2,1)+M(:,1,2).*M(:,3,3).*M(:,4,1).*N(:,2,1)+M(:,1,3).*M(:,3,1).*M(:,4,2).*N(:,2,1)-M(:,1,1).*M(:,3,3).*M(:,4,2).*N(:,2,1)-M(:,1,2).*M(:,3,1).*M(:,4,3).*N(:,2,1)+M(:,1,1).*M(:,3,2).*M(:,4,3).*N(:,2,1)-M(:,1,3).*M(:,2,1).*M(:,4,2).*N(:,3,1)+M(:,1,2).*M(:,2,1).*M(:,4,3).*N(:,3,1)+M(:,1,3).*M(:,2,1).*M(:,3,2).*N(:,4,1)-M(:,1,2).*M(:,2,1).*M(:,3,3).*N(:,4,1)+M(:,2,3).*((-M(:,1,2).*M(:,4,1)+M(:,1,1).*M(:,4,2)).*N(:,3,1)+M(:,3,2).*(M(:,4,1).*N(:,1,1)-M(:,1,1).*N(:,4,1))+M(:,3,1).*(-M(:,4,2).*N(:,1,1)+M(:,1,2).*N(:,4,1)))+M(:,2,2).*((M(:,1,3).*M(:,4,1)-M(:,1,1).*M(:,4,3)).*N(:,3,1)+M(:,3,3).*(-M(:,4,1).*N(:,1,1)+M(:,1,1).*N(:,4,1))+M(:,3,1).*(M(:,4,3).*N(:,1,1)-M(:,1,3).*N(:,4,1))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,4,2) =(M(:,2,1).*M(:,3,3).*M(:,4,2).*N(:,1,2)-M(:,2,1).*M(:,3,2).*M(:,4,3).*N(:,1,2)-M(:,1,3).*M(:,3,2).*M(:,4,1).*N(:,2,2)+M(:,1,2).*M(:,3,3).*M(:,4,1).*N(:,2,2)+M(:,1,3).*M(:,3,1).*M(:,4,2).*N(:,2,2)-M(:,1,1).*M(:,3,3).*M(:,4,2).*N(:,2,2)-M(:,1,2).*M(:,3,1).*M(:,4,3).*N(:,2,2)+M(:,1,1).*M(:,3,2).*M(:,4,3).*N(:,2,2)-M(:,1,3).*M(:,2,1).*M(:,4,2).*N(:,3,2)+M(:,1,2).*M(:,2,1).*M(:,4,3).*N(:,3,2)+M(:,1,3).*M(:,2,1).*M(:,3,2).*N(:,4,2)-M(:,1,2).*M(:,2,1).*M(:,3,3).*N(:,4,2)+M(:,2,3).*((-M(:,1,2).*M(:,4,1)+M(:,1,1).*M(:,4,2)).*N(:,3,2)+M(:,3,2).*(M(:,4,1).*N(:,1,2)-M(:,1,1).*N(:,4,2))+M(:,3,1).*(-M(:,4,2).*N(:,1,2)+M(:,1,2).*N(:,4,2)))+M(:,2,2).*((M(:,1,3).*M(:,4,1)-M(:,1,1).*M(:,4,3)).*N(:,3,2)+M(:,3,3).*(-M(:,4,1).*N(:,1,2)+M(:,1,1).*N(:,4,2))+M(:,3,1).*(M(:,4,3).*N(:,1,2)-M(:,1,3).*N(:,4,2))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,4,3) =(M(:,2,1).*M(:,3,3).*M(:,4,2).*N(:,1,3)-M(:,2,1).*M(:,3,2).*M(:,4,3).*N(:,1,3)-M(:,1,3).*M(:,3,2).*M(:,4,1).*N(:,2,3)+M(:,1,2).*M(:,3,3).*M(:,4,1).*N(:,2,3)+M(:,1,3).*M(:,3,1).*M(:,4,2).*N(:,2,3)-M(:,1,1).*M(:,3,3).*M(:,4,2).*N(:,2,3)-M(:,1,2).*M(:,3,1).*M(:,4,3).*N(:,2,3)+M(:,1,1).*M(:,3,2).*M(:,4,3).*N(:,2,3)-M(:,1,3).*M(:,2,1).*M(:,4,2).*N(:,3,3)+M(:,1,2).*M(:,2,1).*M(:,4,3).*N(:,3,3)+M(:,1,3).*M(:,2,1).*M(:,3,2).*N(:,4,3)-M(:,1,2).*M(:,2,1).*M(:,3,3).*N(:,4,3)+M(:,2,3).*((-M(:,1,2).*M(:,4,1)+M(:,1,1).*M(:,4,2)).*N(:,3,3)+M(:,3,2).*(M(:,4,1).*N(:,1,3)-M(:,1,1).*N(:,4,3))+M(:,3,1).*(-M(:,4,2).*N(:,1,3)+M(:,1,2).*N(:,4,3)))+M(:,2,2).*((M(:,1,3).*M(:,4,1)-M(:,1,1).*M(:,4,3)).*N(:,3,3)+M(:,3,3).*(-M(:,4,1).*N(:,1,3)+M(:,1,1).*N(:,4,3))+M(:,3,1).*(M(:,4,3).*N(:,1,3)-M(:,1,3).*N(:,4,3))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
R(:,4,4) =(M(:,2,1).*M(:,3,3).*M(:,4,2).*N(:,1,4)-M(:,2,1).*M(:,3,2).*M(:,4,3).*N(:,1,4)-M(:,1,3).*M(:,3,2).*M(:,4,1).*N(:,2,4)+M(:,1,2).*M(:,3,3).*M(:,4,1).*N(:,2,4)+M(:,1,3).*M(:,3,1).*M(:,4,2).*N(:,2,4)-M(:,1,1).*M(:,3,3).*M(:,4,2).*N(:,2,4)-M(:,1,2).*M(:,3,1).*M(:,4,3).*N(:,2,4)+M(:,1,1).*M(:,3,2).*M(:,4,3).*N(:,2,4)-M(:,1,3).*M(:,2,1).*M(:,4,2).*N(:,3,4)+M(:,1,2).*M(:,2,1).*M(:,4,3).*N(:,3,4)+M(:,1,3).*M(:,2,1).*M(:,3,2).*N(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*N(:,4,4)+M(:,2,3).*((-M(:,1,2).*M(:,4,1)+M(:,1,1).*M(:,4,2)).*N(:,3,4)+M(:,3,2).*(M(:,4,1).*N(:,1,4)-M(:,1,1).*N(:,4,4))+M(:,3,1).*(-M(:,4,2).*N(:,1,4)+M(:,1,2).*N(:,4,4)))+M(:,2,2).*((M(:,1,3).*M(:,4,1)-M(:,1,1).*M(:,4,3)).*N(:,3,4)+M(:,3,3).*(-M(:,4,1).*N(:,1,4)+M(:,1,1).*N(:,4,4))+M(:,3,1).*(M(:,4,3).*N(:,1,4)-M(:,1,3).*N(:,4,4))))./(M(:,1,2).*M(:,2,4).*M(:,3,3).*M(:,4,1)-M(:,1,2).*M(:,2,3).*M(:,3,4).*M(:,4,1)-M(:,1,1).*M(:,2,4).*M(:,3,3).*M(:,4,2)+M(:,1,1).*M(:,2,3).*M(:,3,4).*M(:,4,2)-M(:,1,2).*M(:,2,4).*M(:,3,1).*M(:,4,3)+M(:,1,1).*M(:,2,4).*M(:,3,2).*M(:,4,3)+M(:,1,2).*M(:,2,1).*M(:,3,4).*M(:,4,3)-M(:,1,1).*M(:,2,2).*M(:,3,4).*M(:,4,3)+M(:,1,4).*(M(:,2,3).*(M(:,3,2).*M(:,4,1)-M(:,3,1).*M(:,4,2))+M(:,2,2).*(-M(:,3,3).*M(:,4,1)+M(:,3,1).*M(:,4,3))+M(:,2,1).*(M(:,3,3).*M(:,4,2)-M(:,3,2).*M(:,4,3)))+M(:,1,2).*M(:,2,3).*M(:,3,1).*M(:,4,4)-M(:,1,1).*M(:,2,3).*M(:,3,2).*M(:,4,4)-M(:,1,2).*M(:,2,1).*M(:,3,3).*M(:,4,4)+M(:,1,1).*M(:,2,2).*M(:,3,3).*M(:,4,4)+M(:,1,3).*(M(:,2,4).*(-M(:,3,2).*M(:,4,1)+M(:,3,1).*M(:,4,2))+M(:,2,2).*(M(:,3,4).*M(:,4,1)-M(:,3,1).*M(:,4,4))+M(:,2,1).*(-M(:,3,4).*M(:,4,2)+M(:,3,2).*M(:,4,4))));
