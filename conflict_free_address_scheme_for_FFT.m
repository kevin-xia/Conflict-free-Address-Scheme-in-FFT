%% The address scheme mainly originates from the following publication:
%% T. Pei-Yun and L. Chung-Yi, ¡°A generalized conflict-free memory addressing scheme for continuous-flow parallel-processing FFT processors with rescheduling,¡± 
%% IEEE Trans. VLSI Syst., vol. 19, no. 12, pp. 2290¨C2302, Dec. 2011.
%% The code is available for free trial, non-commercial research or education purposes, and for non-profit organizations. 
%% If you plan on using the code or parts thereof for commercial purposes or if you intend to re-distribute the code or parts thereof, you must contact the author. 
%% If you are using the code or parts thereof for your scientific work, you must to provide a reference to this website or the above publication.
%% Author: Kaifeng Xia (e-mail: xiakaifeng@ime.ac.cn)

%%------------------------------------------------------
% The test script is for the address distibution method for memoey-based FFT in 2048 points. The decomposition algorithm is radix-2^3 and the PE number is 2.
% It can also be used to test FFTs with other points by changing the address scheme. 
% Decomposition manners of different lengths are: 64--->(2^3)*(2^3); 128--->2*(2^3)*(2^3); 256--->(2^2)*(2^3)*(2^3);
%                                                512--->(2^3)*(2^3)*(2^3);  1024--->2*(2^3)*(2^3)*(2^3);  2048--->(2^2)*(2^3)*(2^3)*(2^3);
%			                         4096--->(2^3)*(2^3)*(2^3)*(2^3);

clear;
clc;
close all;

pnt_inf=[64	128	256	512	1024	2048	4096];

din_num_u_fwd(1,2048)=0;din_num_l_fwd(1,2048)=0;
dout_num_u_fwd(1,2048)=0;dout_num_l_fwd(1,2048)=0;
din_num_u_inv(1,2048)=0;din_num_l_inv(1,2048)=0;
dout_num_u_inv(1,2048)=0;dout_num_l_inv(1,2048)=0;
num_add{1,4}=0;qs_lev(1,4)=0;ms_lev(1,4)=0;
num_bank(1,4096)=0;
num_pnt(1,4096)=0;
pnt_addr0 = 0; pnt_addr1 = 0; pnt_addr2 = 0; pnt_addr3 = 0; 


for num_i=6:6
	len_bit=log2(pnt_inf(num_i));
	len_lev=ceil(len_bit/3);
	rdx_lev1=3*(mod(len_bit,3)==0)+mod(len_bit,3);
	ms=0;
	for num_lev=1:len_lev
		din_num_u_fwd(1,2048)=0;din_num_l_fwd(1,2048)=0;
		dout_num_u_fwd(1,2048)=0;dout_num_l_fwd(1,2048)=0;
		din_num_u_inv(1,2048)=0;din_num_l_inv(1,2048)=0;
		dout_num_u_inv(1,2048)=0;dout_num_l_inv(1,2048)=0;
		if(num_lev==1)
			qs=rdx_lev1;
		else
			qs=3;
		end
		ms=ms+qs;
		for num_cnt=0:pnt_inf(num_i)/2-1
			cnt_bit=de2bi(num_cnt,len_bit-1,'left-msb');
			din_bit_ass_u=[cnt_bit(1:len_bit-1-qs+1) 0 cnt_bit(len_bit-1-qs+2:len_bit-1)];
			din_bit_ass_l=[cnt_bit(1:len_bit-1-qs+1) 1 cnt_bit(len_bit-1-qs+2:len_bit-1)];
			din_bit_u=[din_bit_ass_u(len_bit-ms+1:len_bit) din_bit_ass_u(1:len_bit-ms)];
			din_bit_l=[din_bit_ass_l(len_bit-ms+1:len_bit) din_bit_ass_l(1:len_bit-ms)];
			dout_bit_ass_u=[cnt_bit(1:len_bit-1) 0];
			dout_bit_ass_l=[cnt_bit(1:len_bit-1) 1];
			dout_bit_u=[dout_bit_ass_u(len_bit-ms+1:len_bit) dout_bit_ass_u(1:len_bit-ms)];
			dout_bit_l=[dout_bit_ass_l(len_bit-ms+1:len_bit) dout_bit_ass_l(1:len_bit-ms)];
			din_num_u_fwd(num_cnt+1)=bi2de(din_bit_u,2,'left-msb');
			din_num_l_fwd(num_cnt+1)=bi2de(din_bit_l,2,'left-msb');
			dout_num_u_fwd(num_cnt+1)=bi2de(dout_bit_u,2,'left-msb');
			dout_num_l_fwd(num_cnt+1)=bi2de(dout_bit_l,2,'left-msb');
			din_num_u_inv(num_cnt+1)=bi2de(din_bit_u,2,'right-msb');
			din_num_l_inv(num_cnt+1)=bi2de(din_bit_l,2,'right-msb');
			dout_num_u_inv(num_cnt+1)=bi2de(dout_bit_u,2,'right-msb');
			dout_num_l_inv(num_cnt+1)=bi2de(dout_bit_l,2,'right-msb');
		end
		if(qs==3)
			din_num_u_fwd=reshape(din_num_u_fwd,4,512);din_num_l_fwd=reshape(din_num_l_fwd,4,512);
			dout_num_u_fwd=reshape(dout_num_u_fwd,4,512);dout_num_l_fwd=reshape(dout_num_l_fwd,4,512);
			din_num_u_inv=reshape(din_num_u_inv,4,512);din_num_l_inv=reshape(din_num_l_inv,4,512);
			dout_num_u_inv=reshape(dout_num_u_inv,4,512);dout_num_l_inv=reshape(dout_num_l_inv,4,512);
		elseif(qs==2)
			din_num_u_fwd=reshape(din_num_u_fwd,2,1024);din_num_l_fwd=reshape(din_num_l_fwd,2,1024);
			dout_num_u_fwd=reshape(dout_num_u_fwd,2,1024);dout_num_l_fwd=reshape(dout_num_l_fwd,2,1024);
			din_num_u_inv=reshape(din_num_u_inv,2,1024);din_num_l_inv=reshape(din_num_l_inv,2,1024);
			dout_num_u_inv=reshape(dout_num_u_inv,2,1024);dout_num_l_inv=reshape(dout_num_l_inv,2,1024);
		end
		num.din_up_data_fwd=din_num_u_fwd;num.din_dn_data_fwd=din_num_l_fwd;	%input and output addresses
		num.dout_up_data_fwd=dout_num_u_fwd;num.dout_dn_data_fwd=dout_num_l_fwd;
		num.din_up_data_inv=din_num_u_inv;num.din_dn_data_inv=din_num_l_inv;	%inverse input and output addresses
		num.dout_up_data_inv=dout_num_u_inv;num.dout_dn_data_inv=dout_num_l_inv;
		num_add{num_lev}=num;
		clear	din_num_u_fwd din_num_l_fwd dout_num_u_fwd dout_num_l_fwd din_num_u_inv din_num_l_inv dout_num_u_inv dout_num_l_inv;
		
		qs_lev(num_lev)=qs;ms_lev(num_lev)=ms;
	end
	
	for num_cnt=0:pnt_inf(num_i)-1
		cnt_bit=de2bi(num_cnt,len_bit,'right-msb');			%memory bank index
		bank_sel_1=xor(cnt_bit(1),xor(cnt_bit(2),xor(cnt_bit(3),xor(cnt_bit(5),xor(cnt_bit(7),xor(cnt_bit(9),xor(cnt_bit(10),cnt_bit(11))))))));
		bank_sel_2=xor(cnt_bit(4),xor(cnt_bit(6),cnt_bit(8)));

		num_bank(num_cnt+1)=bank_sel_1*2+bank_sel_2;			%memory bank index of each order
		
		if(num_bank(num_cnt+1)==0)
			pnt_addr0 = pnt_addr0 + 1;
			num_pnt(num_cnt+1) = pnt_addr0;
		elseif(num_bank(num_cnt+1)==1)
			pnt_addr1 = pnt_addr1 + 1;
			num_pnt(num_cnt+1) = pnt_addr1;
		elseif(num_bank(num_cnt+1)==2)
			pnt_addr2 = pnt_addr2 + 1;
			num_pnt(num_cnt+1) = pnt_addr2;
		else
			pnt_addr3 = pnt_addr3 + 1;
			num_pnt(num_cnt+1) = pnt_addr3;
		end
	end



	for num_lev=1:len_lev						       %test the conflictions
		num_tem=num_add{num_lev};
		din_num_u_fwd=num_tem.din_up_data_fwd;din_num_l_fwd=num_tem.din_dn_data_fwd;
		dout_num_u_fwd=num_tem.dout_up_data_fwd;dout_num_l_fwd=num_tem.dout_dn_data_fwd;
		din_num_u_fwd_bank=num_bank(din_num_u_fwd+1);din_num_l_fwd_bank=num_bank(din_num_l_fwd+1);
		dout_num_u_fwd_bank=num_bank(din_num_u_fwd+1);dout_num_l_fwd_bank=num_bank(din_num_l_fwd+1);
		din_num_u_inv=num_tem.din_up_data_inv;din_num_l_inv=num_tem.din_dn_data_inv;
		dout_num_u_inv=num_tem.dout_up_data_inv;dout_num_l_inv=num_tem.dout_dn_data_inv;
		din_num_u_inv_bank=num_bank(din_num_u_inv+1);din_num_l_inv_bank=num_bank(din_num_l_inv+1);
		dout_num_u_inv_bank=num_bank(din_num_u_inv+1);dout_num_l_inv_bank=num_bank(din_num_l_inv+1);
		
 		for num_cnt=0:511
 			bi_num=de2bi(num_cnt,9,'left-msb');
 			
 			if(num_lev==1)
 				din_aaa_bi=[bi_num(1:5) 0 bi_num(6:8) 0 bi_num(9)];dout_aaa_bi=[bi_num(1:5) 0 bi_num(6:9) 0];
 				din_bbb_bi=[bi_num(1:5) 0 bi_num(6:8) 1 bi_num(9)];dout_bbb_bi=[bi_num(1:5) 0 bi_num(6:9) 1];
 				din_ccc_bi=[bi_num(1:5) 1 bi_num(6:8) 0 bi_num(9)];dout_ccc_bi=[bi_num(1:5) 1 bi_num(6:9) 0];
 				din_ddd_bi=[bi_num(1:5) 1 bi_num(6:8) 1 bi_num(9)];dout_ddd_bi=[bi_num(1:5) 1 bi_num(6:9) 1];
 			elseif(num_lev==2)
 				din_aaa_bi=[bi_num(1:2) 0 bi_num(3:7) 0 bi_num(8:9)];dout_aaa_bi=[bi_num(1:2) 0 bi_num(3:9) 0];
 				din_bbb_bi=[bi_num(1:2) 0 bi_num(3:7) 1 bi_num(8:9)];dout_bbb_bi=[bi_num(1:2) 0 bi_num(3:9) 1];
 				din_ccc_bi=[bi_num(1:2) 1 bi_num(3:7) 0 bi_num(8:9)];dout_ccc_bi=[bi_num(1:2) 1 bi_num(3:9) 0];
 				din_ddd_bi=[bi_num(1:2) 1 bi_num(3:7) 1 bi_num(8:9)];dout_ddd_bi=[bi_num(1:2) 1 bi_num(3:9) 1];
 			elseif(num_lev==3)
 				din_aaa_bi=[bi_num(1:7) 0  0 bi_num(8:9)];dout_aaa_bi=[bi_num(1:7) 0 bi_num(8:9) 0];
 				din_bbb_bi=[bi_num(1:7) 0  1 bi_num(8:9)];dout_bbb_bi=[bi_num(1:7) 0 bi_num(8:9) 1];
 				din_ccc_bi=[bi_num(1:7) 1  0 bi_num(8:9)];dout_ccc_bi=[bi_num(1:7) 1 bi_num(8:9) 0];
 				din_ddd_bi=[bi_num(1:7) 1  1 bi_num(8:9)];dout_ddd_bi=[bi_num(1:7) 1 bi_num(8:9) 1];
 			else
 				din_aaa_bi=[bi_num(1:3) 0 bi_num(4:7) 0 bi_num(8:9)];dout_aaa_bi=[bi_num(1:3) 0 bi_num(4:9) 0];
 				din_bbb_bi=[bi_num(1:3) 0 bi_num(4:7) 1 bi_num(8:9)];dout_bbb_bi=[bi_num(1:3) 0 bi_num(4:9) 1];
 				din_ccc_bi=[bi_num(1:3) 1 bi_num(4:7) 0 bi_num(8:9)];dout_ccc_bi=[bi_num(1:3) 1 bi_num(4:9) 0];
 				din_ddd_bi=[bi_num(1:3) 1 bi_num(4:7) 1 bi_num(8:9)];dout_ddd_bi=[bi_num(1:3) 1 bi_num(4:9) 1];
 			end
 			
 			if(num_lev==1)
 				din_aaa_num_fwd(num_cnt+1)=bi2de([din_aaa_bi(10:11) din_aaa_bi(1:9)],'left-msb');dout_aaa_num_fwd(num_cnt+1)=bi2de([dout_aaa_bi(10:11) dout_aaa_bi(1:9)],'left-msb');
 				din_bbb_num_fwd(num_cnt+1)=bi2de([din_bbb_bi(10:11) din_bbb_bi(1:9)],'left-msb');dout_bbb_num_fwd(num_cnt+1)=bi2de([dout_bbb_bi(10:11) dout_bbb_bi(1:9)],'left-msb');
 				din_ccc_num_fwd(num_cnt+1)=bi2de([din_ccc_bi(10:11) din_ccc_bi(1:9)],'left-msb');dout_ccc_num_fwd(num_cnt+1)=bi2de([dout_ccc_bi(10:11) dout_ccc_bi(1:9)],'left-msb');
 				din_ddd_num_fwd(num_cnt+1)=bi2de([din_ddd_bi(10:11) din_ddd_bi(1:9)],'left-msb');dout_ddd_num_fwd(num_cnt+1)=bi2de([dout_ddd_bi(10:11) dout_ddd_bi(1:9)],'left-msb');
 				din_aaa_num_inv(num_cnt+1)=bi2de([din_aaa_bi(10:11) din_aaa_bi(1:9)],'right-msb');dout_aaa_num_inv(num_cnt+1)=bi2de([dout_aaa_bi(10:11) dout_aaa_bi(1:9)],'right-msb');
 				din_bbb_num_inv(num_cnt+1)=bi2de([din_bbb_bi(10:11) din_bbb_bi(1:9)],'right-msb');dout_bbb_num_inv(num_cnt+1)=bi2de([dout_bbb_bi(10:11) dout_bbb_bi(1:9)],'right-msb');
 				din_ccc_num_inv(num_cnt+1)=bi2de([din_ccc_bi(10:11) din_ccc_bi(1:9)],'right-msb');dout_ccc_num_inv(num_cnt+1)=bi2de([dout_ccc_bi(10:11) dout_ccc_bi(1:9)],'right-msb');
 				din_ddd_num_inv(num_cnt+1)=bi2de([din_ddd_bi(10:11) din_ddd_bi(1:9)],'right-msb');dout_ddd_num_inv(num_cnt+1)=bi2de([dout_ddd_bi(10:11) dout_ddd_bi(1:9)],'right-msb');
 			elseif(num_lev==2)
 				din_aaa_num_fwd(num_cnt+1)=bi2de([din_aaa_bi(7:11) din_aaa_bi(1:6)],'left-msb');dout_aaa_num_fwd(num_cnt+1)=bi2de([dout_aaa_bi(7:11) dout_aaa_bi(1:6)],'left-msb');
 				din_bbb_num_fwd(num_cnt+1)=bi2de([din_bbb_bi(7:11) din_bbb_bi(1:6)],'left-msb');dout_bbb_num_fwd(num_cnt+1)=bi2de([dout_bbb_bi(7:11) dout_bbb_bi(1:6)],'left-msb');
 				din_ccc_num_fwd(num_cnt+1)=bi2de([din_ccc_bi(7:11) din_ccc_bi(1:6)],'left-msb');dout_ccc_num_fwd(num_cnt+1)=bi2de([dout_ccc_bi(7:11) dout_ccc_bi(1:6)],'left-msb');
 				din_ddd_num_fwd(num_cnt+1)=bi2de([din_ddd_bi(7:11) din_ddd_bi(1:6)],'left-msb');dout_ddd_num_fwd(num_cnt+1)=bi2de([dout_ddd_bi(7:11) dout_ddd_bi(1:6)],'left-msb');
 				din_aaa_num_inv(num_cnt+1)=bi2de([din_aaa_bi(7:11) din_aaa_bi(1:6)],'right-msb');dout_aaa_num_inv(num_cnt+1)=bi2de([dout_aaa_bi(7:11) dout_aaa_bi(1:6)],'right-msb');
 				din_bbb_num_inv(num_cnt+1)=bi2de([din_bbb_bi(7:11) din_bbb_bi(1:6)],'right-msb');dout_bbb_num_inv(num_cnt+1)=bi2de([dout_bbb_bi(7:11) dout_bbb_bi(1:6)],'right-msb');
 				din_ccc_num_inv(num_cnt+1)=bi2de([din_ccc_bi(7:11) din_ccc_bi(1:6)],'right-msb');dout_ccc_num_inv(num_cnt+1)=bi2de([dout_ccc_bi(7:11) dout_ccc_bi(1:6)],'right-msb');
 				din_ddd_num_inv(num_cnt+1)=bi2de([din_ddd_bi(7:11) din_ddd_bi(1:6)],'right-msb');dout_ddd_num_inv(num_cnt+1)=bi2de([dout_ddd_bi(7:11) dout_ddd_bi(1:6)],'right-msb');
 			elseif(num_lev==3)
 				din_aaa_num_fwd(num_cnt+1)=bi2de([din_aaa_bi(4:11) din_aaa_bi(1:3)],'left-msb');dout_aaa_num_fwd(num_cnt+1)=bi2de([dout_aaa_bi(4:11) dout_aaa_bi(1:3)],'left-msb');
 				din_bbb_num_fwd(num_cnt+1)=bi2de([din_bbb_bi(4:11) din_bbb_bi(1:3)],'left-msb');dout_bbb_num_fwd(num_cnt+1)=bi2de([dout_bbb_bi(4:11) dout_bbb_bi(1:3)],'left-msb');
 				din_ccc_num_fwd(num_cnt+1)=bi2de([din_ccc_bi(4:11) din_ccc_bi(1:3)],'left-msb');dout_ccc_num_fwd(num_cnt+1)=bi2de([dout_ccc_bi(4:11) dout_ccc_bi(1:3)],'left-msb');
 				din_ddd_num_fwd(num_cnt+1)=bi2de([din_ddd_bi(4:11) din_ddd_bi(1:3)],'left-msb');dout_ddd_num_fwd(num_cnt+1)=bi2de([dout_ddd_bi(4:11) dout_ddd_bi(1:3)],'left-msb');
 				din_aaa_num_inv(num_cnt+1)=bi2de([din_aaa_bi(4:11) din_aaa_bi(1:3)],'right-msb');dout_aaa_num_inv(num_cnt+1)=bi2de([dout_aaa_bi(4:11) dout_aaa_bi(1:3)],'right-msb');
 				din_bbb_num_inv(num_cnt+1)=bi2de([din_bbb_bi(4:11) din_bbb_bi(1:3)],'right-msb');dout_bbb_num_inv(num_cnt+1)=bi2de([dout_bbb_bi(4:11) dout_bbb_bi(1:3)],'right-msb');
 				din_ccc_num_inv(num_cnt+1)=bi2de([din_ccc_bi(4:11) din_ccc_bi(1:3)],'right-msb');dout_ccc_num_inv(num_cnt+1)=bi2de([dout_ccc_bi(4:11) dout_ccc_bi(1:3)],'right-msb');
 				din_ddd_num_inv(num_cnt+1)=bi2de([din_ddd_bi(4:11) din_ddd_bi(1:3)],'right-msb');dout_ddd_num_inv(num_cnt+1)=bi2de([dout_ddd_bi(4:11) dout_ddd_bi(1:3)],'right-msb');
 			else
 				din_aaa_num_fwd(num_cnt+1)=bi2de(din_aaa_bi,'left-msb');dout_aaa_num_fwd(num_cnt+1)=bi2de(dout_aaa_bi,'left-msb');
 				din_bbb_num_fwd(num_cnt+1)=bi2de(din_bbb_bi,'left-msb');dout_bbb_num_fwd(num_cnt+1)=bi2de(dout_bbb_bi,'left-msb');
 				din_ccc_num_fwd(num_cnt+1)=bi2de(din_ccc_bi,'left-msb');dout_ccc_num_fwd(num_cnt+1)=bi2de(dout_ccc_bi,'left-msb');
 				din_ddd_num_fwd(num_cnt+1)=bi2de(din_ddd_bi,'left-msb');dout_ddd_num_fwd(num_cnt+1)=bi2de(dout_ddd_bi,'left-msb');
 				din_aaa_num_inv(num_cnt+1)=bi2de(din_aaa_bi,'right-msb');dout_aaa_num_inv(num_cnt+1)=bi2de(dout_aaa_bi,'right-msb');
 				din_bbb_num_inv(num_cnt+1)=bi2de(din_bbb_bi,'right-msb');dout_bbb_num_inv(num_cnt+1)=bi2de(dout_bbb_bi,'right-msb');
 				din_ccc_num_inv(num_cnt+1)=bi2de(din_ccc_bi,'right-msb');dout_ccc_num_inv(num_cnt+1)=bi2de(dout_ccc_bi,'right-msb');
 				din_ddd_num_inv(num_cnt+1)=bi2de(din_ddd_bi,'right-msb');dout_ddd_num_inv(num_cnt+1)=bi2de(dout_ddd_bi,'right-msb');
 			end
 		end
 		
 		if(num_lev==1)
 			len_num=length(din_aaa_num_fwd);
 			din_num_u_fwd=din_num_u_fwd(:,1:len_num);din_num_l_fwd=din_num_l_fwd(:,1:len_num);
 			dout_num_u_fwd=dout_num_u_fwd(:,1:len_num);dout_num_l_fwd=dout_num_l_fwd(:,1:len_num);
 			
 			din_num_u_inv=din_num_u_inv(:,1:len_num);din_num_l_inv=din_num_l_inv(:,1:len_num);
 			dout_num_u_inv=dout_num_u_inv(:,1:len_num);dout_num_l_inv=dout_num_l_inv(:,1:len_num);
 			for num_a=1:len_num
 				[din_aaa_row,din_aaa_col]=find(din_num_u_fwd==din_aaa_num_fwd(num_a));din_aaa_bank_fwd=din_num_u_fwd_bank(din_aaa_row,din_aaa_col);[dout_aaa_row,dout_aaa_col]=find(dout_num_u_fwd==dout_aaa_num_fwd(num_a));dout_aaa_bank_fwd=dout_num_u_fwd_bank(dout_aaa_row,dout_aaa_col);
 				[din_bbb_row,din_bbb_col]=find(din_num_l_fwd==din_bbb_num_fwd(num_a));din_bbb_bank_fwd=din_num_l_fwd_bank(din_bbb_row,din_bbb_col);[dout_bbb_row,dout_bbb_col]=find(dout_num_l_fwd==dout_bbb_num_fwd(num_a));dout_bbb_bank_fwd=dout_num_l_fwd_bank(dout_bbb_row,dout_bbb_col);
 				[din_ccc_row,din_ccc_col]=find(din_num_u_fwd==din_ccc_num_fwd(num_a));din_ccc_bank_fwd=din_num_u_fwd_bank(din_ccc_row,din_ccc_col);[dout_ccc_row,dout_ccc_col]=find(dout_num_u_fwd==dout_ccc_num_fwd(num_a));dout_ccc_bank_fwd=dout_num_u_fwd_bank(dout_ccc_row,dout_ccc_col);
 				[din_ddd_row,din_ddd_col]=find(din_num_l_fwd==din_ddd_num_fwd(num_a));din_ddd_bank_fwd=din_num_l_fwd_bank(din_ddd_row,din_ddd_col);[dout_ddd_row,dout_ddd_col]=find(dout_num_l_fwd==dout_ddd_num_fwd(num_a));dout_ddd_bank_fwd=dout_num_l_fwd_bank(dout_ddd_row,dout_ddd_col);
 				
 				[din_aaa_row,din_aaa_col]=find(din_num_u_inv==din_aaa_num_inv(num_a));din_aaa_bank_inv=din_num_u_inv_bank(din_aaa_row,din_aaa_col);[dout_aaa_row,dout_aaa_col]=find(dout_num_u_inv==dout_aaa_num_inv(num_a));dout_aaa_bank_inv=dout_num_u_inv_bank(dout_aaa_row,dout_aaa_col);
 				[din_bbb_row,din_bbb_col]=find(din_num_l_inv==din_bbb_num_inv(num_a));din_bbb_bank_inv=din_num_l_inv_bank(din_bbb_row,din_bbb_col);[dout_bbb_row,dout_bbb_col]=find(dout_num_l_inv==dout_bbb_num_inv(num_a));dout_bbb_bank_inv=dout_num_l_inv_bank(dout_bbb_row,dout_bbb_col);
 				[din_ccc_row,din_ccc_col]=find(din_num_u_inv==din_ccc_num_inv(num_a));din_ccc_bank_inv=din_num_u_inv_bank(din_ccc_row,din_ccc_col);[dout_ccc_row,dout_ccc_col]=find(dout_num_u_inv==dout_ccc_num_inv(num_a));dout_ccc_bank_inv=dout_num_u_inv_bank(dout_ccc_row,dout_ccc_col);
 				[din_ddd_row,din_ddd_col]=find(din_num_l_inv==din_ddd_num_inv(num_a));din_ddd_bank_inv=din_num_l_inv_bank(din_ddd_row,din_ddd_col);[dout_ddd_row,dout_ddd_col]=find(dout_num_l_inv==dout_ddd_num_inv(num_a));dout_ddd_bank_inv=dout_num_l_inv_bank(dout_ddd_row,dout_ddd_col);
 				if(length(unique([din_aaa_bank_fwd,din_bbb_bank_fwd,din_ccc_bank_fwd,din_ddd_bank_fwd]))~=4)
 					spr_hdl_in=sprintf('the %d forwrd input not meet the rule at level %d.',num_a,num_lev);
 					disp(spr_hdl_in);
 				end
 				
 				if(length(unique([dout_aaa_bank_fwd,dout_bbb_bank_fwd,dout_ccc_bank_fwd,dout_ddd_bank_fwd]))~=4)
 					spr_hdl_out=sprintf('the %d forwrd output not meet the rule at level %d.',num_a,num_lev);
 					disp(spr_hdl_out);
 				end
 				
				if(length(unique([din_aaa_bank_inv,din_bbb_bank_inv,din_ccc_bank_inv,din_ddd_bank_inv]))~=4)
 					spr_hdl_in=sprintf('the %d inverse input not meet the rule at level %d.',num_a,num_lev);
 					disp(spr_hdl_in);
 				end
 				
 				if(length(unique([dout_aaa_bank_inv,dout_bbb_bank_inv,dout_ccc_bank_inv,dout_ddd_bank_inv]))~=4)
 					spr_hdl_out=sprintf('the %d inverse output not meet the rule at level %d.',num_a,num_lev);
 					disp(spr_hdl_out);
 				end
 			end
 		elseif(num_lev==2 || num_lev==3 || num_lev==4)
 			len_num=length(din_aaa_num_fwd);
 			din_num_u_fwd=din_num_u_fwd(:,1:len_num/2);din_num_l_fwd=din_num_l_fwd(:,1:len_num/2);
 			dout_num_u_fwd=dout_num_u_fwd(:,1:len_num/2);dout_num_l_fwd=dout_num_l_fwd(:,1:len_num/2);
 			
 			din_num_u_inv=din_num_u_inv(:,1:len_num/2);din_num_l_inv=din_num_l_inv(:,1:len_num/2);
 			dout_num_u_inv=dout_num_u_inv(:,1:len_num/2);dout_num_l_inv=dout_num_l_inv(:,1:len_num/2);
 			for num_a=1:len_num
 				[din_aaa_row,din_aaa_col]=find(din_num_u_fwd==din_aaa_num_fwd(num_a));din_aaa_bank_fwd=din_num_u_fwd_bank(din_aaa_row,din_aaa_col);[dout_aaa_row,dout_aaa_col]=find(dout_num_u_fwd==dout_aaa_num_fwd(num_a));dout_aaa_bank_fwd=dout_num_u_fwd_bank(dout_aaa_row,dout_aaa_col);
 				[din_bbb_row,din_bbb_col]=find(din_num_l_fwd==din_bbb_num_fwd(num_a));din_bbb_bank_fwd=din_num_l_fwd_bank(din_bbb_row,din_bbb_col);[dout_bbb_row,dout_bbb_col]=find(dout_num_l_fwd==dout_bbb_num_fwd(num_a));dout_bbb_bank_fwd=dout_num_l_fwd_bank(dout_bbb_row,dout_bbb_col);
 				[din_ccc_row,din_ccc_col]=find(din_num_u_fwd==din_ccc_num_fwd(num_a));din_ccc_bank_fwd=din_num_u_fwd_bank(din_ccc_row,din_ccc_col);[dout_ccc_row,dout_ccc_col]=find(dout_num_u_fwd==dout_ccc_num_fwd(num_a));dout_ccc_bank_fwd=dout_num_u_fwd_bank(dout_ccc_row,dout_ccc_col);
 				[din_ddd_row,din_ddd_col]=find(din_num_l_fwd==din_ddd_num_fwd(num_a));din_ddd_bank_fwd=din_num_l_fwd_bank(din_ddd_row,din_ddd_col);[dout_ddd_row,dout_ddd_col]=find(dout_num_l_fwd==dout_ddd_num_fwd(num_a));dout_ddd_bank_fwd=dout_num_l_fwd_bank(dout_ddd_row,dout_ddd_col);
 				
 				[din_aaa_row,din_aaa_col]=find(din_num_u_inv==din_aaa_num_inv(num_a));din_aaa_bank_inv=din_num_u_inv_bank(din_aaa_row,din_aaa_col);[dout_aaa_row,dout_aaa_col]=find(dout_num_u_inv==dout_aaa_num_inv(num_a));dout_aaa_bank_inv=dout_num_u_inv_bank(dout_aaa_row,dout_aaa_col);
 				[din_bbb_row,din_bbb_col]=find(din_num_l_inv==din_bbb_num_inv(num_a));din_bbb_bank_inv=din_num_l_inv_bank(din_bbb_row,din_bbb_col);[dout_bbb_row,dout_bbb_col]=find(dout_num_l_inv==dout_bbb_num_inv(num_a));dout_bbb_bank_inv=dout_num_l_inv_bank(dout_bbb_row,dout_bbb_col);
 				[din_ccc_row,din_ccc_col]=find(din_num_u_inv==din_ccc_num_inv(num_a));din_ccc_bank_inv=din_num_u_inv_bank(din_ccc_row,din_ccc_col);[dout_ccc_row,dout_ccc_col]=find(dout_num_u_inv==dout_ccc_num_inv(num_a));dout_ccc_bank_inv=dout_num_u_inv_bank(dout_ccc_row,dout_ccc_col);
 				[din_ddd_row,din_ddd_col]=find(din_num_l_inv==din_ddd_num_inv(num_a));din_ddd_bank_inv=din_num_l_inv_bank(din_ddd_row,din_ddd_col);[dout_ddd_row,dout_ddd_col]=find(dout_num_l_inv==dout_ddd_num_inv(num_a));dout_ddd_bank_inv=dout_num_l_inv_bank(dout_ddd_row,dout_ddd_col);
 				if(length(unique([din_aaa_bank_fwd,din_bbb_bank_fwd,din_ccc_bank_fwd,din_ddd_bank_fwd]))~=4)
 					spr_hdl_in=sprintf('the %d forwrd input not meet the rule at level %d.',num_a,num_lev);
 					disp(spr_hdl_in);
 				end
 				
 				if(length(unique([dout_aaa_bank_fwd,dout_bbb_bank_fwd,dout_ccc_bank_fwd,dout_ddd_bank_fwd]))~=4)
 					spr_hdl_out=sprintf('the %d forwrd output not meet the rule at level %d.',num_a,num_lev);
 					disp(spr_hdl_out);
 				end
 				
				if(length(unique([din_aaa_bank_inv,din_bbb_bank_inv,din_ccc_bank_inv,din_ddd_bank_inv]))~=4)
 					spr_hdl_in=sprintf('the %d inverse input not meet the rule at level %d.',num_a,num_lev);
 					disp(spr_hdl_in);
 				end
 				
 				if(length(unique([dout_aaa_bank_inv,dout_bbb_bank_inv,dout_ccc_bank_inv,dout_ddd_bank_inv]))~=4)
 					spr_hdl_out=sprintf('the %d inverse output not meet the rule at level %d.',num_a,num_lev);
 					disp(spr_hdl_out);
 				end
			end
 		end
 		
 		clear din_aaa_num_fwd din_bbb_num_fwd din_ccc_num_fwd din_ddd_num_fwd dout_aaa_num_fwd dout_bbb_num_fwd dout_ccc_num_fwd dout_ddd_num_fwd ...
 		din_aaa_num_inv din_bbb_num_inv din_ccc_num_inv din_ddd_num_inv dout_aaa_num_inv dout_bbb_num_inv dout_ccc_num_inv dout_ddd_num_inv;
 		
	end
end