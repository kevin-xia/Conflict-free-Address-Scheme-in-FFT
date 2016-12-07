%% The address scheme mainly originates from the following publications:
%% [1] H. Chen-Fong, C. Yuan, and L. Chen-Yi, ※A generalized mixed-radix algorithm for memory-based FFT processors,§ 
%% IEEE Trans. Circuits Syst.II, vol. 57, no. 1, pp. 26每30, Jan. 2010.
%% [2] J. Chen, J. Hu, S. Lee, and G. E. Sobelman, ※Hardware efficient mixed radix-25/16/9 FFT for LTE systems,§ 
%% IEEE Trans. VLSI Syst., vol. 23, no. 2, pp. 221每229, Feb. 2015.
%% [3] D. Kolba and T. W. Parks, ※A prime factor FFT algorithm using high-speed convolution,§ 
%% IEEE Trans. Acoust., Speech, Signal Process.,vol. 25, no. 4, pp. 281每294, Apr. 1977.
%% [4] C. Temperton, ※A generalized prime factor FFT algorithm for any n = 2p3q5r,§ 
%% SIAM J. on Scientific and Statistical Computing, vol. 13, no. 3,pp. 676每686, Mar. 1992.
%% The code is available for free trial, non-commercial research or education purposes, and for non-profit organizations. 
%% If you plan on using the code or parts thereof for commercial purposes or if you intend to re-distribute the code or parts thereof, you must contact the author. 
%% If you are using the code or parts thereof for your scientific work, you must to provide a reference to this website or the above publications.
%% Author: Kaifeng Xia (e-mail: xiakaifeng@ime.ac.cn)

%%------------------------------------------------------
% The test script is for the address distibution method for memoey-based FFT in 1296 points applied with prime factor algorithm.
% It can also be used to test FFTs with other points by changing the specific address scheme. 
% Decomposition manners of 1296-point is: 1296--->9*9*16(3*3*3*3*4*4);

clear;
clc;
close all;

pnt_num=1296;
num_cnt = 0; num_bank = zeros(1,pnt_num); num_addr = zeros(1,pnt_num);


for lev_1=0:2
	for lev_2=0:2
		for lev_3=0:2
			for lev_4=0:2
				for lev_5=0:3
					for lev_6=0:3
						num_cnt = num_cnt +1;
						%% memory bank and address
						num_bank(num_cnt) = mod((lev_1+lev_2+lev_3+lev_4+lev_5+lev_6),4);
						num_addr(num_cnt) = floor((lev_1*432+lev_2*144+lev_3*48+lev_4*16+lev_5*4+lev_6)/4);
					end
				end
			end
		end
	end
end


for num_lev3=0:143		%level1
	for num_lev2=0:2
		lev_i = 0;
		for num_lev1=0:2
			lev_i = lev_i + 1;
			num_dat(lev_i) = num_lev1*432 + num_lev2*144 + num_lev3;
		end
		bank_sel(1) = num_bank(num_dat(1)+1);
		bank_sel(2) = num_bank(num_dat(2)+1);
		bank_sel(3) = num_bank(num_dat(3)+1);
		
		bank_dect = unique([bank_sel(1),bank_sel(2),bank_sel(3)]);
		
		if(length(bank_dect)~=3)
			spr_hdl=sprintf('the input not meet the rule at level %d' ,1);
			disp(spr_hdl);
		end
	end
end

for num_lev1=0:8		%level2
	for num_lev4=0:15
		for num_lev3=0:2
			lev_i = 0;
			for num_lev2=0:2
				lev_i = lev_i + 1;
				num_dat(lev_i) = num_lev1*144 + num_lev2*48 + num_lev3*16 + num_lev4;
			end
			bank_sel(1) = num_bank(num_dat(1)+1); bank_sel(2) = num_bank(num_dat(2)+1);
			bank_sel(3) = num_bank(num_dat(3)+1);
			bank_dect = unique([bank_sel(1),bank_sel(2),bank_sel(3)]);
			if(length(bank_dect)~=3)
				spr_hdl=sprintf('the input not meet the rule at level %d' ,2);
				disp(spr_hdl);
			end
		end
	end
end

for num_lev1=0:8		%level3
	for num_lev2=0:8
		for num_lev4=0:3
			lev_i = 0;
			for num_lev3=0:3
				lev_i = lev_i + 1;
				num_dat(lev_i) = num_lev1*144 + num_lev2*16 + num_lev3*4+num_lev4;
			end
			bank_sel(1) = num_bank(num_dat(1)+1); bank_sel(2) = num_bank(num_dat(2)+1);
			bank_sel(3) = num_bank(num_dat(3)+1); bank_sel(4) = num_bank(num_dat(4)+1);
			bank_dect = unique([bank_sel(1),bank_sel(2),bank_sel(3),bank_sel(4)]);
			if(length(bank_dect)~=4)
				spr_hdl=sprintf('the input not meet the rule at level %d ',3);
				disp(spr_hdl);
			end
		end
	end
end