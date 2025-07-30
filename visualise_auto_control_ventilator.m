figure()
plot(out.tout,out.adjustedSetPoint.Data,"LineWidth",2)
hold on
plot(out.tout,out.driving_pressure_dyn.Data,"LineWidth",2)

ylabel("Pressure (cmH2O)")
yyaxis right
plot(out.tout,out.flowIndex.Data,".")
ylabel("Flow Index")
legend("PS","P_dyn","FI")
xlabel("Time (s)")