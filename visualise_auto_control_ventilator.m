t_fi =out.tout(find(out.flowIndex.Data>0));
fi =out.flowIndex.Data(find(out.flowIndex.Data>0));
fi_smooth =smooth(fi,5);


figure()
plot(out.tout,out.adjustedSetPoint.Data,"LineWidth",2)
hold on
plot(out.tout,out.driving_pressure_dyn.Data,".");
plot(out.tout(20000:end),smooth(out.driving_pressure_dyn.Data(20000:end),100000),"LineWidth",2,"HandleVisibility","off","Color","r")
ylabel("Pressure (cmH2O)")
yyaxis right
plot(out.tout,out.flowIndex.Data,".")
hold on
plot(t_fi,fi_smooth,"-","LineWidth",2,"HandleVisibility","off")
ylabel("Flow Index")
legend("PS","P_dyn","FI")
xlabel("Time (s)")

%% max effort vs ps
figure()
plot(out.tout,out.adjustedSetPoint.Data,"LineWidth",2)
hold on
plot(out.tout,smooth(out.adjustedSetPoint.Data,100000),"--","LineWidth",2,"HandleVisibility","off")
plot(out.tout,out.maxEffort1.Data,"LineWidth",2)
plot(out.tout,smooth(out.maxEffort1.Data,100000),"--","LineWidth",2,"HandleVisibility","off")
legend("PS","Pmus")
