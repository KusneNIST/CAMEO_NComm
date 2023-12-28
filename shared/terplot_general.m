function [h,hg,htick]=terplot_general(tx, ty)
%FUNCTION [h,hg,htick]=TERPLOT perpares a ternary axis system that is
% needed for the ternaryc function. It returns three handels:
% - h:      to modify the patch created by the fill function;
% - hg:     to change each grid line separately (must probably be
% modified);
% - htick:  to edit the tick labels (probably very inconvinient)
%
% Uli Theune, Geophysics, University of Alberta

% MODIFICATION COPYRIGHT:
% This software was developed by employees of the National Institute of
% Standards and Technology (NIST), an agency of the Federal Government and
% is being made available as a public service. Pursuant to title 17 United
% States Code Section 105, works of NIST employees are not subject to
% copyright protection in the United States.  This software may be subject
% to foreign copyright.  Permission in the United States and in foreign
% countries, to the extent that NIST may hold copyright, to use, copy,
% modify, create derivative works, and distribute this software and its
% documentation without fee is hereby granted on a non-exclusive basis,
% provided that this notice and disclaimer of warranty appears in all
% copies.

% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
% EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
% WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
% WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND
% FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL
% CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR
% FREE.  IN NO EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT
% NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES,
% ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS
% SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR
% OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY OR
% OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT OF
% THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.

% A. Gilad Kusne, NIST, aaron.kusne@nist.gov, Release 8/01/2020
% If using this work for a publication, please cite:
% Kusne, A. G., et al. "On-the-fly Closed-loop Autonomous Materials
% Discovery via Bayesian Active Learning." arXiv preprint arXiv:2006.06141
% (2020).
%%
%tx = [0 .6 .3 0]';
ty = [0 0 .3*sqrt(3) 0]';
number=11;

%h=fill([0 1 0.5 0],[0 0 0.866 0],'w','linewidth',2);
h=fill(tx,ty,'w','linewidth',2);
%set(h,'facecolor',[0.7 0.7 0.7],'edgecolor','w')
%set(gcf,'color',[0 0 0.3])
d1=cos(pi/3);
d2=sin(pi/3);
l=linspace(0,1,number);
gray = [170,170,170]/255;
hold on
for i=2:length(l)-1
    x1 = linspace(l(i)*d1, 1-l(i)*d1,100);
    y1 = linspace(l(i)*d2, l(i)*d2, 100);
    x2 = linspace(l(i), l(i)+(1-l(i))*d1, 100);
    y2 = linspace(0, (1-l(i))*d2, 100);
    x3 = linspace((1-l(i))*d1, 1-l(i), 100);
    y3 = linspace((1-l(i))*d2, 0, 100);
    in1 = inpolygon(x1,y1,tx,ty);
    in2 = inpolygon(x2,y2,tx,ty);
    in3 = inpolygon(x3,y3,tx,ty);
   plot(x1(in1),y1(in1),'color',gray,'linewidth',1);
   plot(x2(in2),y2(in2),'color',gray,'linewidth',1);
   plot(x3(in3),y3(in3),'color',gray,'linewidth',1);
end
for i=1:(length(tx)-1)
    line([tx(i) tx(i+1)],[ty(i) ty(i+1)],'color','k','linewidth',2);
end
hold off
axis image
axis off
% Make x-tick labels
% for i=1:number
%     htick(i,1)=text(l(i),-0.025,num2str(l(i)));
%     htick(i,3)=text(1-l(i)*cos(pi/3)+0.025,l(i)*sin(pi/3)+0.025,num2str(l(i)));
%     htick(i,2)=text(0.5-l(i)*cos(pi/3)-0.06,sin(pi/3)*(1-l(i)),num2str(l(i)));
% end
htick = [];
%set(gcf,'WindowButtonDownFcn','InitTerExpl');
global hx;



