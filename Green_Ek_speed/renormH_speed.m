function [Hs,Hb]=renormH_speed(E,H,T,Ni)
	%iteration process
% 	Hs=H+T/(E-H)*T';
% 	Hb=H+T/(E-H)*T'+T'/(E-H)*T;
% 	Tt=T/(E-H)*T;
% 	Ttp=T'/(E-H)*T';
% 	for i=1:Ni-1
% 		Hst=Hs+Tt/(E-Hb)*Ttp;
% 		Hbt=Hb+Tt/(E-Hb)*Ttp+Ttp/(E-Hb)*Tt;
% 		Ttt=Tt/(E-Hb)*Tt;
% 		Ttpt=Ttp/(E-Hb)*Ttp;
% 		Hs=Hst;
% 		Hb=Hbt;
% 		Tt=Ttt;
% 		Ttp=Ttpt;
%   end
    tmp_1 = T/(E-H);
    tmp_2 = T'/(E-H);
    tmp_3 = tmp_1*T';
	Hs=H+tmp_3;
	Hb=H+tmp_3+tmp_2*T;
	Tt=tmp_1*T;
	Ttp=tmp_2*T';
    for i=1:Ni-1
        tmp_4 = Tt/(E-Hb);
        tmp_5 = Ttp/(E-Hb);
        tmp_6 = tmp_4*Ttp;
		Hst=Hs+tmp_6;
		Hbt=Hb+tmp_6+tmp_5*Tt;
		Ttt=tmp_4*Tt;
		Ttpt=tmp_5*Ttp;
		Hs=Hst;
		Hb=Hbt;
		Tt=Ttt;
		Ttp=Ttpt;
    end
end