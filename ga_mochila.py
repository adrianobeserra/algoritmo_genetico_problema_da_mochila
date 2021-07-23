from enum import Enum
import copy
import random

tamanho_populacao = 4
porcentagem_mutacao = 0.10
massa_maxima_mochila = 25

prob_acumuladas_tipo1 = {0: 0.25, 1: 0.5, 2: 0.75, 3: 1}
prob_acumuladas_tipo2 = {0: 0.333, 1: 0.667, 2: 1}
prob_acumuladas_tipo3 = {0: 0.167, 1: 0.333,
                         2: 0.5, 3: 0.667, 4: 0.834, 5: 1}

probabilidades_genes = (0.333, 0.666, 1.0)

class AlgoritmoGenetico:
    def __init__(self) -> None:
        super().__init__()
        self.__populacao = []
        self.__melhor = None
        self.__contador_geracao = 1
        self.__iteracoes_sem_melhor_solucao = 0

    def run(self):
        ag.inicializar_populacao()
        while self.__contador_geracao < 150 and self.__iteracoes_sem_melhor_solucao < 10:
            avaliacoes = ag.avaliar_populacao()
            paiA, paiB = ag.selecionar_pais(avaliacoes)
            filho1, filho2 = ag.crossover(paiA, paiB)
            ag.mutacao(filho1, filho2)
            ag.imprimir_populacao()
        print("Execucao encerrada. Melhor solução: ", self.__melhor[0], ". Melhor Cromossomo: ", self.__melhor[1])
        print("Quantidade de geracoes: ", self.__contador_geracao)

    def inicializar_populacao(self):
        print("Populacao Inicial")
        while(len(self.__populacao) < tamanho_populacao):
            quantidade_gene1 = self.__determinar_quantidade_gene(prob_acumuladas_tipo1)
            gene1 = Gene(quantidade_gene1, Item.ITEM_1)

            quantidade_gene2 = self.__determinar_quantidade_gene(prob_acumuladas_tipo2)
            gene2 = Gene(quantidade_gene2, Item.ITEM_2)

            quantidade_gene3 = self.__determinar_quantidade_gene(prob_acumuladas_tipo3)
            gene3 = Gene(quantidade_gene3, Item.ITEM_3)

            novo_cromossomo = Cromossomo(gene1, gene2, gene3)

            if (CromossomoUtils.cromossomo_eh_valido(novo_cromossomo)):
                self.__populacao.append(novo_cromossomo)
                print(novo_cromossomo)
            else:
                print(novo_cromossomo,  " <- invalido")

    def avaliar_populacao(self):
        avaliacoes = []
        for individuo in self.__populacao:
            avaliacoes.append((CromossomoUtils.obter_funcao_objetivo(individuo), individuo))
        melhor_atual_populacao = max(avaliacoes, key=lambda item: item[0])
        if (self.__melhor == None):
            self.__melhor = melhor_atual_populacao
            print("Melhor solução encontrada: ", melhor_atual_populacao[0], melhor_atual_populacao[1])
        elif melhor_atual_populacao[0] > self.__melhor[0]:
            self.__melhor = melhor_atual_populacao
            self.__iteracoes_sem_melhor_solucao = 0
            print("Melhor solução encontrada: ", melhor_atual_populacao[0], melhor_atual_populacao[1])
        else:
            self.__iteracoes_sem_melhor_solucao = self.__iteracoes_sem_melhor_solucao + 1
            print("Solução inferior encontrada: ", melhor_atual_populacao[0], melhor_atual_populacao[1])
            
            
        return avaliacoes

    def selecionar_pais(self, avaliacoes):
        avaliacoes_ordenadas = sorted(avaliacoes, key=lambda item: item[0], reverse=True)
        indice_meio = len(avaliacoes_ordenadas) // 2

        melhores_avaliacoes = avaliacoes_ordenadas[:indice_meio]
        melhores_avaliacoes = self.__obter_probabilidade_acumulada(melhores_avaliacoes)
        paiA = self.__sortear_item(melhores_avaliacoes)

        piores_avaliacoes = avaliacoes_ordenadas[indice_meio:]
        piores_avaliacoes = self.__obter_probabilidade_acumulada(piores_avaliacoes)
        paiB = self.__sortear_item(piores_avaliacoes)

        print('Pais: ' + '[' + ', '.join([str(gene.quantidade_item) for gene in paiA.genes]) +
              ']' + '[' + ', '.join([str(gene.quantidade_item) for gene in paiB.genes]) + ']')

        return paiA, paiB

    def crossover(self, paiA, paiB):
        novos_pais_sao_validos = False
        filho1 = None
        filho2 = None

        while novos_pais_sao_validos == False:
            gene_a_trocar1 = self.__sortear_gene()
            gene_a_trocar2 = self.__sortear_gene()

            filho1 = copy.deepcopy(paiA)
            filho2 = copy.deepcopy(paiB)

            gene_aux = filho1.genes[gene_a_trocar1]
            filho1.genes[gene_a_trocar1] = filho2.genes[gene_a_trocar1]
            filho2.genes[gene_a_trocar1] = gene_aux

            if (gene_a_trocar2 != gene_a_trocar1):
                gene_aux = filho1.genes[gene_a_trocar2]
                filho1.genes[gene_a_trocar2] = filho2.genes[gene_a_trocar2]
                filho2.genes[gene_a_trocar2] = gene_aux

            if CromossomoUtils.cromossomo_eh_valido(filho1) and CromossomoUtils.cromossomo_eh_valido(filho2) and filho1 != paiA and filho2 != paiB:
                novos_pais_sao_validos = True

                print('Filhos: ' + '[' + ', '.join([str(gene.quantidade_item) for gene in filho1.genes]
                                                   ) + ']' + '[' + ', '.join([str(gene.quantidade_item) for gene in filho2.genes]) + ']')

                self.__populacao.remove(paiA)
                self.__populacao.remove(paiB)
                self.__populacao.append(filho1)
                self.__populacao.append(filho2)

        return filho1, filho2

    def mutacao(self, individuo1, individuo2):
        numero_aleatorio = self.__obter_numero_aleatorio()
        if (numero_aleatorio <= porcentagem_mutacao):
            numero_aleatorio = self.__obter_numero_aleatorio()
            if (numero_aleatorio <= 0.5):
                individuo_pre_mutacao = individuo1
            else:
                individuo_pre_mutacao = individuo2
            gene_a_trocar = self.__sortear_gene()

            novo_gene_eh_valido = False
            while novo_gene_eh_valido == False:
                novo_gene = None
                if (gene_a_trocar == 0):
                    nova_quantidade_gene = self.__determinar_quantidade_gene(prob_acumuladas_tipo1)
                    novo_gene = Gene(nova_quantidade_gene, Item.ITEM_1)
                elif (gene_a_trocar == 1):
                    nova_quantidade_gene = self.__determinar_quantidade_gene(prob_acumuladas_tipo2)
                    novo_gene = Gene(nova_quantidade_gene, Item.ITEM_2)
                else:
                    nova_quantidade_gene = self.__determinar_quantidade_gene(prob_acumuladas_tipo3)
                    novo_gene = Gene(nova_quantidade_gene, Item.ITEM_3)
                
                individuo_pos_mutacao = copy.deepcopy(individuo_pre_mutacao)
                individuo_pos_mutacao.genes[gene_a_trocar] = novo_gene
                if (CromossomoUtils.cromossomo_eh_valido(individuo_pos_mutacao) and individuo_pos_mutacao != individuo_pre_mutacao):
                    novo_gene_eh_valido = True
                    self.__populacao[self.__populacao.index(individuo_pre_mutacao)] = individuo_pos_mutacao
                    print('Filho com mutação sofrida: ', individuo_pre_mutacao, ' -> ', individuo_pos_mutacao)
    
    def imprimir_populacao(self):
        self.__contador_geracao = self.__contador_geracao + 1
        print("\n##### Populacao ", self.__contador_geracao)
        for cromossomo in self.__populacao:
            print(cromossomo)
        

    def __obter_numero_aleatorio(self):
        return random.randint(1, 99)/100

    def __obter_probabilidade_acumulada(self, avaliacoes):
        avaliacoes_com_probabilidade_acumulada = []
        total_melhores_avaliacoes = sum(avaliacao[0] for avaliacao in avaliacoes)
        probabilidade_acumulada = 0
        for avaliacao in avaliacoes:
            probabilidade_acumulada += (avaliacao[0] / total_melhores_avaliacoes)
            avaliacao = avaliacao + (probabilidade_acumulada,)
            avaliacoes_com_probabilidade_acumulada.append(avaliacao)
        return avaliacoes_com_probabilidade_acumulada

    def __sortear_item(self, probabilidades_acumuladas):
        numero_aleatorio = self.__obter_numero_aleatorio()
        for cont, tupla in enumerate(probabilidades_acumuladas):
            if (cont < len(probabilidades_acumuladas) - 1):
                indice_prob_acumulada_atual = cont
                indice_prob_acumulada_proxima = cont + 1

                if (numero_aleatorio <= probabilidades_acumuladas[indice_prob_acumulada_atual][2]):
                    return probabilidades_acumuladas[indice_prob_acumulada_atual][1]
                elif (numero_aleatorio > probabilidades_acumuladas[indice_prob_acumulada_atual][2] and numero_aleatorio < probabilidades_acumuladas[indice_prob_acumulada_proxima][2]):
                    return probabilidades_acumuladas[indice_prob_acumulada_proxima][1]

        prob_acumulada_ultima = len(probabilidades_acumuladas) - 1
        return probabilidades_acumuladas[prob_acumulada_ultima][1]

    def __sortear_gene(self):
        numero_aleatorio = self.__obter_numero_aleatorio()
        for cont, tupla in enumerate(probabilidades_genes):
            if (cont < len(probabilidades_genes) - 1):
                indice_prob_acumulada_atual = cont
                indice_prob_acumulada_proxima = cont

                if (numero_aleatorio <= probabilidades_genes[indice_prob_acumulada_atual]):
                    return indice_prob_acumulada_atual
                elif (numero_aleatorio > probabilidades_genes[indice_prob_acumulada_atual] and numero_aleatorio < probabilidades_genes[indice_prob_acumulada_proxima]):
                    return indice_prob_acumulada_proxima

        indice_prob_acumulada_ultima = len(probabilidades_genes) - 1
        return indice_prob_acumulada_ultima

    def __determinar_quantidade_gene(self, probabilidades_acumuladas):
        numero_aleatorio = self.__obter_numero_aleatorio()
        for cont, total in enumerate(probabilidades_acumuladas):
            if (cont < len(probabilidades_acumuladas) - 1):
                prob_acumulada_atual = total
                prob_acumulada_proxima = total + 1

                if (numero_aleatorio <= probabilidades_acumuladas[prob_acumulada_atual]):
                    return prob_acumulada_atual
                elif (numero_aleatorio > probabilidades_acumuladas[prob_acumulada_atual] and numero_aleatorio < probabilidades_acumuladas[prob_acumulada_proxima]):
                    return prob_acumulada_proxima

        prob_acumulada_ultima = len(probabilidades_acumuladas) - 1
        return probabilidades_acumuladas[prob_acumulada_ultima]


class Cromossomo:
    def __init__(self, gene1, gene2, gene3):
        self.__genes = [gene1, gene2, gene3]

    @property
    def genes(self):
        return self.__genes

    @genes.setter
    def genes(self, genes):
        self.__genes = genes

    def __str__(self):
        return '[' + ', '.join([str(gene.quantidade_item) for gene in self.__genes]) + ']'

    def __eq__(self, other):
        valores_outros_genes = [gene.quantidade_item for gene in other.genes]
        valores_genes_atual = [gene.quantidade_item for gene in self.genes]
        if valores_outros_genes == valores_genes_atual:
            return True
        else:
            return False


class CromossomoUtils:
    def __init__(self) -> None:
        super().__init__()

    @staticmethod
    def cromossomo_eh_valido(cromossomo):
        if (CromossomoUtils.__obter_massa_cromossomo(cromossomo) >= massa_maxima_mochila):
            return False
        else:
            return True

    @staticmethod
    def __obter_massa_cromossomo(cromossomo):
        return sum(gene.quantidade_item * gene.item.massa for gene in cromossomo.genes)

    @staticmethod
    def obter_funcao_objetivo(cromossomo):
        return sum(gene.quantidade_item * gene.item.valor for gene in cromossomo.genes)


class Gene:
    def __init__(self, quantidade_item, item):
        self.__quantidade_item = quantidade_item
        self.__item = item

    @property
    def quantidade_item(self):
        return self.__quantidade_item

    @quantidade_item.setter
    def quantidade_item(self, quantidade_item):
        self.__quantidade_item = quantidade_item

    @property
    def item(self):
        return self.__item

    @item.setter
    def item(self, item):
        self.__item = item


class Item(Enum):
    ITEM_1 = 1, 3, 40
    ITEM_2 = 2, 5, 100
    ITEM_3 = 3, 2, 50

    def __new__(cls, *args, **kwds):
        value = len(cls.__members__) + 1
        obj = object.__new__(cls)
        obj._value_ = value
        return obj

    def __init__(self, tipo, massa, valor):
        self.tipo = tipo
        self.massa = massa
        self.valor = valor

if __name__ == "__main__":
    ag = AlgoritmoGenetico()
    ag.run()