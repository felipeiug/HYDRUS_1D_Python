# Projeto de Conversão Fortran → Python

Este projeto tem como objetivo **converter um código legado em Fortran** para **Python**, garantindo maior portabilidade, facilidade de manutenção e integração com bibliotecas modernas como **NumPy** e **Pandas**.

---

## 📂 Estrutura do Projeto

- `H1D_Fortran/`  
  Contém o código original em **Fortran**, que servirá de referência para a conversão.

- `src/`  
  Implementação em **Python** equivalente ao código Fortran.  
  Aqui estarão os módulos reescritos e testados.

- `README.md`  
  Este arquivo com instruções e informações do projeto.

---

## 🎯 Objetivo

- Preservar a **lógica numérica e científica** do código original em Fortran.
- Traduzir as rotinas para Python, aproveitando recursos como:
  - **NumPy** para operações matriciais e vetoriais.
  - **Pandas** para manipulação de dados em tabelas.
  - **Matplotlib** para visualização dos resultados.
- Garantir compatibilidade de resultados entre as duas versões.

---

## ⚙️ Pré-requisitos

Certifique-se de ter instalado:

- [Python 3.10+](https://www.python.org/downloads/)
- [NumPy](https://numpy.org/)  
- [Pandas](https://pandas.pydata.org/)  
- [Matplotlib](https://matplotlib.org/)  

Certifique-se de utilizar uma virtualenv para organização das dependências executando:
```bash
python -m virtualenv venv
venv/Scripts/activate
```

Instale as dependências executando:

```bash
pip install -r requirements.txt
```

---

## ▶️ Como Executar

1. **Rodar a versão Python**:

```bash
python src/main.py
```

2. (Opcional) Comparar saídas entre Fortran e Python.

---

## 🛠️ Metodologia de Conversão

1. **Análise** do código em `H1D_Fortran/`.
2. **Mapeamento** das funções e estruturas Fortran para equivalentes Python.
3. **Reescrita** progressiva dos módulos.
4. **Validação** comparando resultados numéricos.
5. **Documentação** do processo e das funções implementadas para futuras manutenções.

---

## 📌 Observações

- O código Fortran original é mantido apenas como **referência histórica** e para validação.  
- A versão Python será o código **principal e atualizado**.  
- Sempre que possível, funções específicas do `MSFLIB` ou outras bibliotecas legadas serão substituídas por equivalentes Python.

---

## 👨‍💻 Autores

Projeto de conversão e manutenção realizado por **Felipeiug**.
