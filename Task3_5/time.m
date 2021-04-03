N_x = [10 20 40 80 160];
t_GS_f = [6.9700000000000014E-004  8.1449999999999995E-003  0.10383600000000000       1.4015979999999999    19.594902999999999]            
t_GS_c = [2.1900000000000001E-004  1.6269999999999995E-003  2.1097999999999999E-002   0.31472999999999995   5.1984430000000001   ]
t_J_c  = [4.4600000000000022E-004  2.9979999999999998E-003  4.4144000000000003E-002   0.69909900000000003   12.787506000000000] 
t_J_f  = [2.4380000000000001E-003  2.1913999999999999E-002  0.19643300000000000       2.6312829999999998    41.715816000000004  ]

figure(1)
axis square
set (gca,'position',[0.15,0.15,0.8,0.8] );
F1 = loglog(N_x,t_GS_f,'-D','DisplayName','Gauss Seidel (Tol $1e-12$)','linewidth',2)
hold on
F2 = loglog(N_x,t_GS_c,'-*','DisplayName','Gauss Seidel (Tol $h_x^2$)','linewidth',2)
F3 = loglog(N_x,t_J_f,'--p','DisplayName','Jacobi (Tol $1e-12$)','linewidth',2)
F4 = loglog(N_x,t_J_c,'.-.','DisplayName','Jacobi (Tol $h_x^2$)','linewidth',2)
F5 = loglog(N_x(2:3),2e-7*N_x(2:3).^4,'DisplayName','slope 4','linewidth',2)
hold off
legend([F1,F2,F3,F4,F5],'Interpreter','latex','Location','northwest')
xlabel("$\log_{10}$N",'Interpreter','latex')
ylabel("$\log10$ CPU time",'Interpreter','latex')
set(gca, 'FontSize', 20)
grid on
print Figure_Neumann_BC_cpu_time.eps -depsc2 -r600
