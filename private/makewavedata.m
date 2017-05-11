function eta=makewavedata(eamp,ffreqs,wns,dirs,layout,datatypes,depth,fs,ndat)

%makewavedata: make pseudo random sea elevation data for a specified layout of probes
%
%[eta] = makerandomsea(eamp,wns,dirs,xx,yy)
%
%outputs:
%eta		is output surface elevation
%
%inputs
%eamp		matrix of the component amplitudes [i,j] with i'th frequency component and j'th directional component
%ffreqs	column vector of component radian frequencies. These must correspond to the frequency bins in input eamp
%wns 		column vector of component wavenumbers in m-1. These must correspond to the frequency bins in input eamp
%dirs		row vector of component directions in rads. These must correspond to the directional bins in input eamp
%layout	matrix containing x,y,z coordinates of instrument location with format:
%				[x1 x2 x3...]
%				[y1 y2 y3...]
%				[z1 z2 z3...])
%data type matrix containing sensor types for DIWASP transfer functions: must have same number of columns as layout
%
%fs		sampling frequency
%ndat		length of data
%
%amplitude matrix must be [m by n] where m is number of wavenumber 
%components in wns and n is number of directional components in dirs
%

%Copyright (C) 2002 Coastal Oceanography Group, CWR, UWA, Perth


disp('generating instrument data');

randphase=8*atan(1.0)*rand(size(eamp));
nprobes=size(layout,2);
dt=1/fs;
t=[0:dt:dt*ndat];

freqmat=ffreqs*ones(size(dirs));

%get transfer functions for instrument type
for i=1:nprobes
    trf(:,:,i)=feval(datatypes{i},ffreqs,dirs,wns,layout(3,i),depth);
end

%make random data
for i=1:ndat
   for j=1:nprobes
      kxwt=layout(1,j)*wns*cos(dirs) + layout(2,j)*wns*sin(dirs) - freqmat*t(i)+randphase;
      eta(i,j)=sum(sum( real( eamp.*trf(:,:,j).*exp(-1i*kxwt) ) ));
   end
end


