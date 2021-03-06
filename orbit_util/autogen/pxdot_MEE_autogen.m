function pxdot = pxdot_MEE_autogen(L,T,f,g,h,k,m,mu,p,px1,px2,px3,px4,px5,px6)
%pxdot_MEE_autogen
%    PXDOT = pxdot_MEE_autogen(L,T,F,G,H,K,M,MU,P,PX1,PX2,PX3,PX4,PX5,PX6)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    29-Apr-2022 00:52:48

t2 = cos(L);
t3 = sin(L);
t4 = mu.*p;
t5 = h.^2;
t6 = k.^2;
t17 = 1.0./m;
t18 = 1.0./mu;
t19 = 1.0./p.^2;
t7 = t2.^2;
t8 = t3.^2;
t9 = f.*t2;
t10 = g.*t2;
t11 = h.*t2;
t12 = k.*t2;
t13 = f.*t3;
t14 = g.*t3;
t15 = h.*t3;
t16 = k.*t3;
t27 = p.*t18;
t28 = sqrt(t4);
t29 = t5+t6+1.0;
t20 = t9.*2.0;
t21 = t10.*2.0;
t22 = t13.*2.0;
t23 = t14.*2.0;
t24 = -t13;
t26 = -t15;
t30 = t11+t16;
t31 = t9+t14+1.0;
t33 = sqrt(t27);
t25 = -t22;
t32 = t31+1.0;
t34 = t10+t24;
t35 = conj(t33);
t36 = t12+t26;
t38 = t20+t23+2.0;
t39 = t31.^2;
t42 = 1.0./t31;
t37 = t21+t25;
t40 = t2.*t32;
t41 = t3.*t32;
t43 = 1.0./t39;
t44 = 1.0./t35;
t47 = 1.0./t38;
t49 = px3.*t2.*t35;
t50 = px2.*t3.*t35;
t52 = p.*px1.*t35.*t42.*2.0;
t57 = px6.*t35.*t36.*t42;
t60 = f.*px3.*t35.*t36.*t42;
t61 = g.*px2.*t35.*t36.*t42;
t45 = f+t40;
t46 = g+t41;
t48 = t47.^2;
t51 = -t50;
t58 = px4.*t2.*t29.*t35.*t47;
t59 = px5.*t3.*t29.*t35.*t47;
t62 = -t57;
t63 = -t60;
t53 = t49+t51;
t64 = px2.*t35.*t42.*t45;
t65 = px3.*t35.*t42.*t46;
t70 = t58+t59+t61+t62+t63;
t54 = abs(t53);
t55 = sign(t53);
t66 = t52+t64+t65;
t71 = abs(t70);
t72 = sign(t70);
t56 = t54.^2;
t67 = abs(t66);
t68 = sign(t66);
t73 = t71.^2;
t69 = t67.^2;
t74 = t56+t69+t73;
t75 = 1.0./sqrt(t74);
mt1 = [(T.*t17.*t75.*(t71.*t72.*(px6.*t18.*t36.*t42.*t44.*(-1.0./2.0)-(f.*px3.*t18.*t36.*t42.*t44)./2.0+(g.*px2.*t18.*t36.*t42.*t44)./2.0+(px4.*t2.*t18.*t29.*t44.*t47)./2.0+(px5.*t3.*t18.*t29.*t44.*t47)./2.0).*2.0+t67.*t68.*(px1.*t35.*t42.*2.0+px1.*t27.*t42.*t44+(px2.*t18.*t42.*t44.*t45)./2.0+(px3.*t18.*t42.*t44.*t46)./2.0).*2.0-t54.*t55.*((px2.*t3.*t18.*t44)./2.0-(px3.*t2.*t18.*t44)./2.0).*2.0))./2.0+1.0./p.^3.*px6.*t28.*t39.*2.0-(mu.*px6.*t19.*t39)./(t28.*2.0)];
mt2 = [T.*t17.*t75.*(t71.*t72.*(px3.*t35.*t36.*t42-px6.*t2.*t35.*t36.*t43+px4.*t7.*t29.*t35.*t48.*2.0+px2.*t10.*t35.*t36.*t43-px3.*t9.*t35.*t36.*t43+px5.*t2.*t3.*t29.*t35.*t48.*2.0).*2.0+t67.*t68.*(-t3.*t42.*t49+t43.*t46.*t49-px2.*t35.*t42.*(t7+1.0)+p.*px1.*t2.*t35.*t43.*2.0+px2.*t2.*t35.*t43.*t45).*2.0).*(-1.0./2.0)-px6.*t2.*t19.*t28.*t31.*2.0];
mt3 = [(T.*t17.*t75.*(t71.*t72.*(px2.*t35.*t36.*t42+px6.*t3.*t35.*t36.*t43-px5.*t8.*t29.*t35.*t48.*2.0-px2.*t14.*t35.*t36.*t43+px3.*t13.*t35.*t36.*t43-px4.*t2.*t3.*t29.*t35.*t48.*2.0).*2.0-t67.*t68.*(t2.*t42.*t51+t43.*t45.*t50-px3.*t35.*t42.*(t8+1.0)+p.*px1.*t3.*t35.*t43.*2.0+px3.*t3.*t35.*t43.*t46).*2.0))./2.0-px6.*t3.*t19.*t28.*t31.*2.0;T.*t17.*t71.*t72.*t75.*(px6.*t3.*t35.*t42-px2.*t14.*t35.*t42+px3.*t13.*t35.*t42+px4.*t11.*t35.*t47.*2.0+px5.*t15.*t35.*t47.*2.0);T.*t17.*t71.*t72.*t75.*(-px6.*t2.*t35.*t42+px2.*t10.*t35.*t42-px3.*t9.*t35.*t42+px4.*t12.*t35.*t47.*2.0+px5.*t16.*t35.*t47.*2.0)];
mt4 = [T.*t17.*t75.*(t71.*t72.*(px6.*t30.*t35.*t42+f.*px3.*t30.*t35.*t42-g.*px2.*t30.*t35.*t42-px4.*t3.*t29.*t35.*t47+px5.*t2.*t29.*t35.*t47+px6.*t34.*t35.*t36.*t43+f.*px3.*t34.*t35.*t36.*t43-g.*px2.*t34.*t35.*t36.*t43-px4.*t2.*t29.*t35.*t37.*t48-px5.*t3.*t29.*t35.*t37.*t48).*-2.0+t54.*t55.*(px2.*t2.*t35+px3.*t3.*t35).*2.0+t67.*t68.*(px2.*t35.*t42.*(t41-t2.*t34)-px3.*t35.*t42.*(t40+t3.*t34)+p.*px1.*t34.*t35.*t43.*2.0+px2.*t34.*t35.*t43.*t45+px3.*t34.*t35.*t43.*t46).*2.0).*(-1.0./2.0)-px6.*t19.*t28.*t31.*t34.*2.0];
pxdot = [mt1;mt2;mt3;mt4];
