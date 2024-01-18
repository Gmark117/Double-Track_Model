function plot_data_raw(tyre_data, cut_start, cut_end, name)
    
    FZ = tyre_data.FZ;
    IA = tyre_data.IA;
    SA = tyre_data.SA;
    SL = tyre_data.SL;
    P = tyre_data.P;
    TSTC = tyre_data.TSTC;
    TSTI = tyre_data.TSTI;
    TSTO = tyre_data.TSTO;

    figure('Name', name)
    tiledlayout(6,1)
    
    ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
    plot(-FZ)
    hold on
    plot([cut_start cut_start],y_range,'--r')
    plot([cut_end cut_end],y_range,'--r')
    title('Vertical force')
    xlabel('Samples [-]')
    ylabel('[N]')
    
    ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
    plot(IA)
    hold on
    plot([cut_start cut_start],y_range,'--r')
    plot([cut_end cut_end],y_range,'--r')
    title('Camber angle')
    xlabel('Samples [-]')
    ylabel('[deg]')
    
    ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
    plot(SA)
    hold on
    plot([cut_start cut_start],y_range,'--r')
    plot([cut_end cut_end],y_range,'--r')
    title('Side slip')
    xlabel('Samples [-]')
    ylabel('[deg]')
    
    ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
    plot(SL)
    hold on
    plot([cut_start cut_start],y_range,'--r')
    plot([cut_end cut_end],y_range,'--r')
    title('Longitudinal slip')
    xlabel('Samples [-]')
    ylabel('[-]')
    
    ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
    plot(P)
    hold on
    plot([cut_start cut_start],y_range,'--r')
    plot([cut_end cut_end],y_range,'--r')
    title('Tyre pressure')
    xlabel('Samples [-]')
    ylabel('[psi]')
    
    ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
    plot(TSTC,'DisplayName','Center')
    hold on
    plot(TSTI,'DisplayName','Internal')
    plot(TSTO,'DisplayName','Outboard')
    hold on
    plot([cut_start cut_start],y_range,'--r')
    plot([cut_end cut_end],y_range,'--r')
    title('Tyre temperatures')
    xlabel('Samples [-]')
    ylabel('[degC]')
    
    linkaxes(ax_list,'x')

end