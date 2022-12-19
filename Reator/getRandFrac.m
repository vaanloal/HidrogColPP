function [frac] = getRandFrac()

    % AA(1)   DIA(2) TEA(3) PAA(4)  DEA(5) NEA(6)  PIA(7) DPIAs(8)
    frac_article = [0.487772819 0.046512298 0.028606858 0.164908937 0.035763384 0.138086399 0.098349306 0];
    mu_articles = [36.02742646	4.977268992	5.814219925	25.972811	5.921397934	15.71761666	9.869994036	0.265];
    std_articles = [12.39355057	0.344659155	7.177229778	9.889254176	2.424230593	4.563799252	0.939610102	0.374766594];
    min_articles = [10.88 4.621806977 0.6 15 3.03 5	8.56 0];
    max_articles = [56 5.31 14 54.18 10	27 10.67 0.53];
    mediana_articles = [37	5	2.842659775	24	5	15.72363636	10.12498807	0.265];

    pd1 = makedist('Triangular', 'A', min_articles(1), 'B', mu_articles(1), 'C', max_articles(1));
    pd2 = makedist('Triangular', 'A', min_articles(2), 'B', mu_articles(2), 'C', max_articles(2));
    pd3 = makedist('Triangular', 'A', min_articles(3), 'B', mu_articles(3), 'C', max_articles(3));
    pd4 = makedist('Triangular', 'A', min_articles(4), 'B', mu_articles(4), 'C', max_articles(4));
    pd5 = makedist('Triangular', 'A', min_articles(5), 'B', mu_articles(5), 'C', max_articles(5));
    pd6 = makedist('Triangular', 'A', min_articles(6), 'B', mu_articles(6), 'C', max_articles(6));
    pd7 = makedist('Triangular', 'A', min_articles(7), 'B', mu_articles(7), 'C', max_articles(7));
    pd8 = makedist('Triangular', 'A', min_articles(8), 'B', mu_articles(8), 'C', max_articles(8));

    frac = zeros(1, length(frac_article));
    frac(1) = random(pd1, 1, 1);
    frac(2) = random(pd2, 1, 1);
    frac(3) = random(pd3, 1, 1);
    frac(4) = random(pd4, 1, 1);
    frac(5) = random(pd5, 1, 1);
    frac(6) = random(pd6, 1, 1);
    frac(7) = random(pd7, 1, 1);
    frac(8) = random(pd8, 1, 1);

    frac = frac / (sum(frac));
end
