# all-paths

Código de exercício para explorar os comportamentos de paths e walks em uma rede viária urbana.

## Notes

- data-raw/ tem os dados de input necessários para o projeto (Bacia do Rio Tamanduateí e zonas OD da região metropolitana de SP)

- Usar o networkx com nx.parallel para realizar os calculos nos grafos

- Usar o matplotlib para visualização, com a skill /dataviz do Claude

## Steps

1. Pegar apenas a porção paulistana da bacia do rio Tamanduateí (dentro da cidade de SP).

2. Montar a área de estudo com base nas zonas OD que fazem interseção com essa porção paulistana (1. ) da bacia do Tamanduateí

3. Baixar a rede viária com o osmnx (vias tipo "drive")

4. Selecionar dois nós aleatórios com base em um parâmetro de shortest path (default = 5)

5. Registrar e quantificar os paths entre os dois nós, com o comprimento mínimo sendo o shortest path e o máximo parametrizado (default = mínimo + 10)

6. Repetir o 5., só que com walks no lugar de paths

7. Montar um scatterplot com as duas distribuições

8. Montar um mapa para cada um (walks x paths)