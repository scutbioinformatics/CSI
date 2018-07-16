function purity = purityMeasure( gnd, res )

a = gnd;
b = res;
M = crosstab(a,b); % you can use also use "confusionmat"
nc = sum(M,1);
mc = max(M,[],1);
purity = sum(mc(nc>0))/sum(nc);
end

