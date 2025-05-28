import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
import pingouin as pg

dados = pd.read_csv('psd_analysis.csv') 

print("\nEstatísticas descritivas:")
print(dados.describe())

plt.figure(figsize=(10, 6))
sns.boxplot(data=dados, palette="Set2")
sns.stripplot(data=dados, color="black", alpha=0.5, jitter=True)
plt.ylabel(r'Potência na Banda $\beta$ (mV²)')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('plots/box_plot.png', bbox_inches="tight")
plt.close()

plt.figure(figsize=(10, 6))
for i in range(len(dados)):
    plt.plot(dados.columns, dados.iloc[i], 'o-', color='gray', alpha=0.4)
plt.plot(dados.columns, dados.median(), 'o-', color='red', linewidth=2, label='Mediana')
plt.ylabel('Valores')
plt.legend()
plt.xticks(rotation=45)
plt.tight_layout()
plt.close()

# 3. Teste de Friedman (alternativa não-paramétrica para ANOVA pareada)
print("\n=== Teste de Friedman ===")
friedman_result = stats.friedmanchisquare(dados['normal'], dados['DP'], dados['DP + ECP'])
print(f"Estatística do teste: {friedman_result.statistic:.4f}")
print(f"Valor-p: {friedman_result.pvalue:.4f}")

if friedman_result.pvalue < 0.05:
    print("\nExistem diferenças significativas entre pelo menos um par de condições (p < 0.05).")
    
    # Testes post-hoc de Wilcoxon pareado com correção de Bonferroni
    print("\n=== Testes Post-hoc de Wilcoxon pareado com correção de Bonferroni ===")

    # Transformar os dados para o formato longo
    dados_long = dados.melt(var_name='Condicao', value_name='Valor')
    dados_long['ID'] = np.tile(np.arange(len(dados)), 3)  # Adiciona identificador por sujeito

    # Executa o teste de Wilcoxon pareado com correção de Bonferroni
    posthoc_results = pg.pairwise_tests(
        data=dados_long,
        dv='Valor',
        within='Condicao',
        subject='ID',
        parametric=False,
        padjust='bonf'
    )

    # Filtra para manter apenas comparações entre condições
    posthoc_results = posthoc_results[posthoc_results['Contrast'] == 'Condicao']

    # Exibe os resultados
    print(posthoc_results[['A', 'B', 'p-unc', 'p-corr']].rename(
        columns={'p-unc': 'p-valor', 'p-corr': 'p-valor (corrigido)'}))

    # Interpretação
    print("\nInterpretação:")
    for _, row in posthoc_results.iterrows():
        sig = "***" if row['p-corr'] < 0.001 else "**" if row['p-corr'] < 0.01 else "*" if row['p-corr'] < 0.05 else "n.s."
        print(f"{row['A']} vs {row['B']}: p = {row['p-corr']:.4f} {sig}")
