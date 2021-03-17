N_x = [20 40 80 160 320];
t_GS_f = [1.3540000000000002E-003  2.1089000000000000E-002 0.29172799999999999    4.4526100000000000    67.885990000000007]
T_GS_f = [1000 3842 14710 56168 213930]             
t_GS_c = [2.4499999999999999E-004  2.8750000000000000E-003 4.7811999999999993E-002  0.87985999999999998    14.377520000000001 ]
T_GS_c = [126 562 2479 10829 46948]

t_J_c  = [3.9299999999999969E-004  6.8539999999999998E-003 0.10710099999999999    1.7280570000000000     29.476431000000002 ] 
T_J_c  = [286 1259 5488 23755 102214]
t_J_f  = [ 2.7309999999999999E-003  4.0849999999999997E-002  0.56385099999999999    8.5593899999999987   134.20466800000000]
T_J_f  = [2033 7820 29952 114399 436201]
figure(1)
axis square
set (gca,'position',[0.15,0.15,0.8,0.8] );
F1 = loglog(N_x,t_GS_f,'-D','DisplayName','Gauss Seidel (Tol $h_x^2$)','linewidth',2)
hold on
F2 = loglog(N_x,t_GS_c,'-*','DisplayName','Gauss Seidel (Tol $1e-12$)','linewidth',2)
F3 = loglog(N_x,t_J_f,'--p','DisplayName','Jacobi (Tol $1e-12$)','linewidth',2)
F4 = loglog(N_x,t_J_c,'.-.','DisplayName','Jacobi (Tol $1e-12$)','linewidth',2)
F5 = loglog(N_x(2:3),5e-8*N_x(2:3).^4,'DisplayName','slope 4','linewidth',2)
hold off
legend([F1,F2,F3,F4,F5],'Interpreter','latex','Location','northwest')
xlabel("$\log_{10}$N",'Interpreter','latex')
ylabel("$\log10$ CPU time",'Interpreter','latex')
set(gca, 'FontSize', 20)
grid on
print Figure_Dichlet_BC_cpu_time.eps -depsc2 -r600

figure(2)
axis square
set (gca,'position',[0.15,0.15,0.8,0.8] );
F1 = loglog(N_x,T_GS_f,'-D','DisplayName','Gauss Seidel (Tol $h_x^2$)','linewidth',2)
hold on
F2 = loglog(N_x,T_GS_c,'-*','DisplayName','Gauss Seidel (Tol $1e-12$)','linewidth',2)
F3 = loglog(N_x,T_J_f,'--p','DisplayName','Jacobi (Tol $1e-12$)','linewidth',2)
F4 = loglog(N_x,T_J_c,'.-.','DisplayName','Jacobi (Tol $1e-12$)','linewidth',2)
F5 = loglog(N_x(2:3),6*N_x(2:3).^2,'DisplayName','slope 2','linewidth',2)
hold off
legend([F1,F2,F3,F4,F5],'Interpreter','latex','Location','northwest')
xlabel("$\log_{10}$N",'Interpreter','latex')
ylabel("$\log10$ Iterative number",'Interpreter','latex')
set(gca, 'FontSize', 20)
grid on
print Figure_Dichlet_BC_iter_num.eps -depsc2 -r600