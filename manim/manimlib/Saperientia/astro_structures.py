from manimlib import *
import numpy as np
from manimlib.Saperientia.Variables import *

def angulo_diedro_de_vetores_cartesianos(A, B, V):
    """
    O que faz
    Calcula o ângulo diedro entre dois planos definidos pelo vértice V e os vetores posição de A e B na esfera.
    Esse ângulo corresponde ao ângulo interno do triângulo esférico no vértice V, útil em trigonometria esférica.

    Como usar
        C = angulo_diedro_de_vetores_cartesianos(A, B, V)
    onde A, B e V são arrays NumPy com coordenadas cartesianas. A função normaliza os vetores, aplica a lei dos quatro elementos esférica e retorna o ângulo em graus.

    Parâmetros
        A np.ndarray Vetor posição do primeiro ponto da esfera.
        B np.ndarray Vetor posição do segundo ponto da esfera.
        V np.ndarray Vetor posição do vértice onde o ângulo será calculado.

    Retorno
        float Ângulo esférico interno no vértice V em graus.

    Observações
        Os vetores são normalizados internamente e o cosseno é protegido com clamp para estabilidade numérica.
    """
    # Normalizar os vetores posição (caso não estejam normalizados)
    A = A / np.linalg.norm(A)
    B = B / np.linalg.norm(B)
    V = V / np.linalg.norm(V)
    
    # Calcular os ângulos esféricos (distâncias angulares sobre a esfera entre os pontos)
    a = np.arccos(np.clip(np.dot(B, V), -1.0, 1.0))  # Ângulo entre B e V 
    b = np.arccos(np.clip(np.dot(A, V), -1.0, 1.0))  # Ângulo entre A e V
    c = np.arccos(np.clip(np.dot(A, B), -1.0, 1.0))  # Ângulo entre A e B
    
    # Aplicar a Lei dos Quatro Elementos Esférica
    cos_C = (np.cos(c) - np.cos(a) * np.cos(b)) / (np.sin(a) * np.sin(b))
    C = np.arccos(np.clip(cos_C, -1.0, 1.0))
    
    # Converter para graus
    return np.degrees(C)

VMOBJECT = True

#ELEMENTOS BÁSICOS
class EsferaCeleste(Sphere):
    """
    O que faz
    Cria uma esfera 3D que representa a esfera celeste, servindo como base para posicionar astros, arcos e círculos.

    Como usar
        esfera = EsferaCeleste(raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=BLUE, resolucao=(80, 80), opacidade=0.3)

    Parâmetros
        raio (float, opcional): Raio da esfera. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        cor (Color, opcional): Cor de preenchimento da esfera. Padrão é BLUE.
        resolucao (tuple, opcional): Resolução da malha na forma linhas, colunas. Padrão é 80, 80.
        opacidade (float, opcional): Opacidade do preenchimento da esfera. Padrão é 0.3.
        **kwargs: Argumentos adicionais herdados de Sphere.

    Observações
        A cor e a opacidade são aplicadas ao preenchimento e funcionam bem como pano de fundo para visualizações astronômicas.
    """
    def __init__(self, raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=BLUE, resolucao=(80, 80), opacidade=0.3, **kwargs):
        # Chama o construtor da classe pai (Sphere)
        super().__init__(radius=raio, resolution=resolucao, **kwargs)
        
        # Armazena os atributos personalizados
        self.cor = cor
        self.opacidade = opacidade

        # Configura as propriedades adicionais
        self.set_opacity(self.opacidade)

        # Configura a cor de preenchimento da esfera com padrão quadriculado (checkerboard)
        self.set_color(cor)


        # Configura o contorno da esfera (stroke) com cor, largura e opacidade especificados

class SuperficieObservador(Disk3D):
    """
    O que faz
    Cria um disco horizontal no plano XY que representa a superfície do observador, como um piso ou horizonte.

    Como usar
        piso = SuperficieObservador(raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=GREEN, opacidade=1)

    Parâmetros
        raio (float, opcional): Raio do disco. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        cor (Color, opcional): Cor do disco. Padrão é GREEN.
        resolucao (tuple, opcional): Resolução da malha do disco. Padrão é 100, 100.
        opacidade (float, opcional): Opacidade do disco. Padrão é 1.
        **kwargs: Argumentos adicionais herdados de Disk3D.
    """
    def __init__(self, raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=GREEN,resolucao=(100,100),opacidade=1, **kwargs):
        self.raio = raio
        super().__init__(resolution=resolucao, radius=raio,**kwargs)
        self.set_color(cor)
        self.set_opacity(opacidade)

class SuperficieObservador(Surface):
    """
    O que faz
    Gera uma superfície circular plana no plano XY por parametrização polar, útil para representar a base de referência do observador.

    Como usar
        piso = SuperficieObservador(raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=GREEN, resolucao=(50, 50), opacidade=1)

    Parâmetros
        raio (float, opcional): Raio do disco paramétrico. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        cor (Color, opcional): Cor do preenchimento. Padrão é GREEN.
        resolucao (tuple, opcional): Resolução da superfície. Padrão é 50, 50.
        opacidade (float, opcional): Opacidade do preenchimento. Padrão é 1.
        **kwargs: Argumentos adicionais herdados de Surface.

    Observações
        A parametrização define u como ângulo e v como raio adimensional entre zero e um.
    """
    def __init__(self, raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=GREEN,resolucao=(50,50),opacidade=1, **kwargs):
        self.raio = raio
        super().__init__(
            u_range=[0, TAU],
            v_range=[0, 1],
            color=cor,
            resolution=resolucao,
            **kwargs
        ) 
    def uv_func(self, u, v):
        return np.array([
            self.raio * np.cos(u) * v,
            self.raio * np.sin(u) * v,
            0
        ]) 
  
class PontoAstro(Sphere):
    """
    O que faz
    Cria um ponto na superfície da esfera celeste a partir de altura e azimute, considerando o norte como o eixo Y positivo.

    Como usar
        astro = PontoAstro(altura=30, azimute=120, tamanho=0.05, cor=WHITE, raio=RENDER_CELESTIAL_SPHERE_RADIUS)

    Parâmetros
        altura (float, opcional): Latitude esférica em graus entre menos noventa e mais noventa. Padrão é 0.
        azimute (float, opcional): Longitude esférica em graus entre zero e trezentos e sessenta a partir do norte para o leste. Padrão é 0.
        tamanho (float, opcional): Raio visual do ponto. Padrão é 0.8 vezes ASTRO_SIZE.
        cor (Color, opcional): Cor do ponto. Padrão é WHITE.
        raio (float, opcional): Raio da esfera de referência para posicionamento. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        center (np.ndarray, opcional): Vetor para deslocar o centro da esfera. Padrão é None.
        **kwargs: Argumentos adicionais herdados de Sphere.

    Observações
        Internamente converte para coordenadas cartesianas e posiciona a pequena esfera no local correto da casca esférica.
    """
    def __init__(self, altura, azimute, tamanho=0.8*ASTRO_SIZE, cor=WHITE, raio=RENDER_CELESTIAL_SPHERE_RADIUS, center=None, **kwargs):
        
        # Armazena os atributos do ponto
        self.altura = altura  # Equivalente à latitude esférica
        self.azimute = azimute  # Equivalente à longitude esférica
        self.cor = cor
        self.raio = raio
        self.tamanho = tamanho

        # Ajusta a longitude para alinhar com o sistema de coordenadas da esfera
        longitude = -self.azimute + 90
            
        # Converte a altura e a longitude para radianos
        theta = np.radians(90 - self.altura)  # Ângulo polar (90° - altura para alinhar com a convenção esférica)
        phi = np.radians(longitude)  # Ângulo azimutal
        
        # Converte as coordenadas esféricas para coordenadas cartesianas
        x = self.raio * np.sin(theta) * np.cos(phi)  # Coordenada X
        y = self.raio * np.sin(theta) * np.sin(phi)  # Coordenada Y
        z = self.raio * np.cos(theta)  # Coordenada Z
        
        # Chama o construtor da classe pai (Dot3D) para criar o ponto na posição calculada
        super().__init__(radius=self.tamanho, color=self.cor, **kwargs)
        
        if center is not None:
            self.move_to(np.array([x, y, z])+center)
        else:
            self.move_to(np.array([x, y, z]))

class P(PontoAstro):
    """
    O que faz
    Fornece um atalho de escrita para PontoAstro com os mesmos parâmetros.

    Como usar
        p = P(altura=20, azimute=45)

    Parâmetros
        Iguais aos de PontoAstro.
    """
    def __init__(self, altura, azimute, tamanho=0.8*ASTRO_SIZE, cor=WHITE, raio=RENDER_CELESTIAL_SPHERE_RADIUS, center=None, **kwargs):
        super().__init__(altura, azimute, tamanho=tamanho, cor=cor, raio=raio,center=center, **kwargs)

class PontoAstroEquatorial(PontoAstro):
    """
    O que faz
    Posiciona um ponto com coordenadas equatoriais de declinação e ascensão reta, convertido ao referencial local por latitude e tempo sideral local.

    Como usar
        estrela = PontoAstroEquatorial(declinacao=15, ascencao_reta_graus=120, latitude=-23, TSL_graus=30)

    Parâmetros
        declinacao (float): Declinação em graus. Padrão é 0.
        ascencao_reta_graus (float): Ascensão reta em graus. Padrão é 0.
        latitude (float): Latitude geográfica do observador em graus. Padrão é 0.
        TSL_graus (float, opcional): Tempo sideral local em graus. Padrão é 0.
        tamanho (float, opcional): Tamanho do ponto. Padrão é 0.8 vezes ASTRO_SIZE.
        cor (Color, opcional): Cor do ponto. Padrão é WHITE.
        raio (float, opcional): Raio da esfera de referência. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        **kwargs: Argumentos adicionais herdados de PontoAstro.

    Observações
        Aplica rotação final no eixo X para alinhar a esfera equatorial com a latitude do observador.
    """
    def __init__(self, declinacao, ascencao_reta_graus, latitude, TSL_graus=0,  tamanho=0.8*ASTRO_SIZE, cor=WHITE, raio=RENDER_CELESTIAL_SPHERE_RADIUS, **kwargs):
        
        # Armazena os atributos do ponto
        self.declinacao = declinacao  # Equivalente à latitude esférica
        self.ascencao_reta = ascencao_reta_graus  # Equivalente à longitude esférica
        self.cor = cor
        self.raio = raio
        self.tamanho = tamanho


        # Converte a altura e a longitude para radianos
         # Ângulo polar (90° - altura para alinhar com a convenção esférica)
        phi = -self.ascencao_reta + TSL_graus + 180  # Ângulo azimutal
        
        # Converte as coordenadas esféricas para coordenadas cartesianas
        
        # Chama o construtor da classe pai (Dot3D) para criar o ponto na posição calculada
        super().__init__(self.declinacao,phi,tamanho=self.tamanho, cor=self.cor, raio=raio,**kwargs)
        
        self.rotate_about_origin((latitude - 90) * DEGREES, X_AXIS)

class PontoVernal(PontoAstroEquatorial):
    """
    O que faz
    Marca o ponto do equinócio vernal para o observador dado, com declinação zero e ascensão reta zero, ajustado por latitude e tempo sideral local.

    Como usar
        gamma = PontoVernal(latitude=-23, TSL_graus=15)

    Parâmetros
        latitude (float, opcional): Latitude do observador em graus. Padrão é 0.
        TSL_graus (float, opcional): Tempo sideral local em graus. Padrão é 0.
        tamanho (float, opcional): Tamanho visual do ponto. Padrão é 0.8 vezes ASTRO_SIZE.
        cor (Color, opcional): Cor do marcador. Padrão é PINK.
        raio (float, opcional): Raio da esfera de referência. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        **kwargs: Argumentos adicionais herdados de PontoAstroEquatorial.
    """
    def __init__(self, latitude, TSL_graus=0,  tamanho=0.8*ASTRO_SIZE, cor=PINK, raio=RENDER_CELESTIAL_SPHERE_RADIUS, **kwargs):
        super().__init__(0,0,latitude=latitude,TSL_graus=TSL_graus,tamanho=tamanho, cor=cor, raio=raio,**kwargs)

##VMOBJECTS

class GrandeArco(VMobject):
    """
    O que faz
    Traça o arco de grande círculo entre dois pontos de superfície usando interpolação esférica.

    Como usar
        arco = GrandeArco(ponto1, ponto2, cor=WHITE, num_pontos=100, center=ORIGIN, espessura=LINE_SIZE)

    Parâmetros
        ponto1 (VMobject): Primeiro ponto sobre a esfera. Padrão é None.
        ponto2 (VMobject): Segundo ponto sobre a esfera. Padrão é None.
        cor (Color, opcional): Cor do traço do arco. Padrão é WHITE.
        num_pontos (int, opcional): Número de amostras do arco. Padrão é 50.
        center (np.ndarray, opcional): Centro da esfera caso deslocada. Padrão é ORIGIN.
        espessura (float, opcional): Largura do traço. Padrão é LINE_SIZE.
    """
    def __init__(self, ponto1, ponto2, cor=WHITE, num_pontos=50, center = None, espessura=LINE_SIZE):
        # Inicializa o objeto como um VGroup (grupo de vetores gráficos do Manim)
        super().__init__()
        if center is None:
            center = ORIGIN
        else:
            center = center
        # Obtém as coordenadas cartesianas dos pontos na esfera
        inicio = ponto1.get_center() - center
        fim = ponto2.get_center() - center
        
        # Verifica se os pontos são iguais; se forem, não há arco a ser desenhado
        if np.allclose(inicio, fim):
            return
        
        # Calcula o ângulo entre os vetores que representam os pontos na esfera
        angulo = np.arccos(np.dot(inicio, fim) / (np.linalg.norm(inicio) * np.linalg.norm(fim)))
        
        # Lista para armazenar os pontos do arco interpolado
        pontos_arco = []
        
        # Gera pontos intermediários ao longo do arco usando interpolação esférica (Slerp)
        for t in range(num_pontos + 1):
            ponto_slerp = (
                np.sin((1 - t / num_pontos) * angulo) * inicio +
                np.sin((t / num_pontos) * angulo) * fim
             ) / np.sin(angulo) + center # Ajusta para o raio da esfera
            pontos_arco.append(ponto_slerp)
        
        # Cria um objeto VMobject para representar o arco
        self.set_points_as_corners(pontos_arco)  # Define os pontos do arco
        self.set_color(cor)
        self.set_stroke(width=espessura)
        self.apply_depth_test()
        
class GrandeCirculo(ParametricCurve):
    """
    O que faz
    Desenha um grande círculo em uma esfera a partir de um vetor normal ao plano do círculo.

    Como usar
        gc = GrandeCirculo(vetor_normal=np.array([0, 1, 0]), raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=YELLOW, espessura=LINE_SIZE)

    Parâmetros
        vetor_normal (np.ndarray, opcional): Vetor unitário normal ao plano do círculo. Padrão é vetor unitário em Z.
        raio (float, opcional): Raio da esfera usada como referência. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        cor (Color, opcional): Cor do traço. Padrão é YELLOW.
        espessura (float, opcional): Largura do traço. Padrão é LINE_SIZE.

    Observações
        O círculo é gerado no plano XY e depois rotacionado para alinhar com a normal desejada.
    """
    def __init__(self, vetor_normal, raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=YELLOW, espessura=LINE_SIZE):
        vetor_normal = vetor_normal / np.linalg.norm(vetor_normal)
        raio_corrigido = raio 
        # Círculo base no plano XY
        def circulo_base(t):
            return raio_corrigido * np.array([np.cos(t), np.sin(t), 0])
        
        super().__init__(lambda t: circulo_base(t),
            t_range=[0, TAU,0.2] )

        # Normaliza o vetor normal


        # Cria o grande círculo
        self.set_color(cor).set_stroke(width=espessura)
        

        # Alinha o círculo com o vetor normal usando .rotate()
        # Calcula o ângulo de rotação entre o vetor [0, 0, 1] e o vetor normal
        angle = np.arccos(np.dot(np.array([0, 0, 1]), vetor_normal))
        axis_of_rotation = np.cross(np.array([0, 0, 1]), vetor_normal)

        # Rotaciona o círculo para alinhar com o vetor normal
        self.rotate(angle, axis_of_rotation)
        self.apply_depth_test()

        # Adiciona a função à instância da classe

class Paralelo(VMobject):
    """
    O que faz
    Desenha um círculo paralelo de latitude fixa sobre a esfera, com possibilidade de rotação segundo a latitude do observador.

    Como usar
        paralelo = Paralelo(altitude=30, latitude=0, raio=RENDER_CELESTIAL_SPHERE_RADIUS)

    Parâmetros
        altitude (float, opcional): Latitude esférica do círculo em graus. Padrão é 0.
        latitude (float, opcional): Latitude do observador em graus para inclinar o plano do círculo. Padrão é 0.
        raio (float, opcional): Raio da esfera de referência. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        num_pontos (int, opcional): Número de pontos para amostragem. Padrão é 400.
        espessura (float, opcional): Largura do traço. Padrão é LINE_SIZE.
        cor (Color, opcional): Cor do traço. Padrão é BLUE_C.
    """
    def __init__(self, altitude, latitude=0, raio=RENDER_CELESTIAL_SPHERE_RADIUS, num_pontos=400, espessura=LINE_SIZE, cor=BLUE_C):

        super().__init__()

        # Converter latitude de graus para radianos
        latitude_rad = np.radians(altitude)

        # Gerar pontos ao longo do círculo de latitude variando a longitude (theta)
        pontos_círculo = []
        for theta in np.linspace(0, 2 * np.pi, num_pontos):
            # A posição no círculo em coordenadas cartesianas
            x = raio * np.cos(latitude_rad) * np.cos(theta)
            y = raio * np.cos(latitude_rad) * np.sin(theta)
            z = raio * np.sin(latitude_rad)

            # Adiciona o ponto ao círculo
            pontos_círculo.append(np.array([x, y, z]))

        # Criar o círculo como um VMobject
        self.set_width(espessura).set_points_as_corners(pontos_círculo).set_color(cor)
        
        # Adicionar o círculo ao VGroup
        if latitude != 0:
            self.rotate_about_origin((latitude - 90) * DEGREES, X_AXIS)
        self.apply_depth_test()

class ArcoParalelo(Paralelo):
    """
    O que faz
    Desenha um arco pertencente a um paralelo de latitude, entre azimutes inicial e final.

    Como usar
        arco = ArcoParalelo(altitude=20, angulo_inicial=0, angulo_final=90, raio=RENDER_CELESTIAL_SPHERE_RADIUS)

    Parâmetros
        altitude (float, opcional): Latitude do paralelo em graus. Padrão é 0.
        angulo_inicial (float, opcional): Azimute inicial em graus. Padrão é 0.
        angulo_final (float, opcional): Azimute final em graus. Padrão é 0.
        raio (float, opcional): Raio da esfera de referência. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        num_pontos (int, opcional): Número de pontos do arco. Padrão é 50.
        espessura (float, opcional): Largura do traço. Padrão é LINE_SIZE.
        cor (Color, opcional): Cor do arco. Padrão é WHITE.
    """
    def __init__(self, altitude, angulo_inicial, angulo_final, raio=RENDER_CELESTIAL_SPHERE_RADIUS, num_pontos=50,espessura=LINE_SIZE, cor=WHITE):
        super().__init__(altitude=altitude, raio=raio, num_pontos=num_pontos, cor=cor)
        angulo_inicio = 90 -angulo_inicial
        angulo_fim = 90-angulo_final
        # Converter graus para radianos
        latitude_rad = np.radians(altitude)
        angulo_inicio_rad = np.radians(angulo_inicio)
        angulo_fim_rad = np.radians(angulo_fim)

        # Gerar pontos apenas no intervalo especificado
        pontos_arco = []
        for theta in np.linspace(angulo_inicio_rad, angulo_fim_rad, num_pontos):
            x = raio * np.cos(latitude_rad) * np.cos(theta)
            y = raio * np.cos(latitude_rad) * np.sin(theta)
            z = raio * np.sin(latitude_rad)
            pontos_arco.append(np.array([x, y, z]))

        # Criar o arco como um VMobject
        self.set_points_as_corners(pontos_arco).set_color(cor).set_width(espessura)
        
class MeridianoLocal(Arc):
    """
    O que faz
    Desenha um semi círculo máximo que representa o meridiano local do observador, passando por zênite e nadir.

    Como usar
        mer = MeridianoLocal(raio=RENDER_CELESTIAL_SPHERE_RADIUS, espessura=LINE_SIZE, cor=ORANGE, opacidade=1)

    Parâmetros
        raio (float, opcional): Raio do semi círculo. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        espessura (float, opcional): Largura do traço. Padrão é LINE_SIZE.
        cor (Color, opcional): Cor do traço. Padrão é ORANGE.
        opacidade (float, opcional): Opacidade do traço. Padrão é 1.
        **kwargs: Argumentos adicionais herdados de Arc.
    """
    def __init__(self, raio=RENDER_CELESTIAL_SPHERE_RADIUS, espessura=LINE_SIZE, cor=ORANGE, opacidade=1, **kwargs):

        super().__init__(
            arc_center=ORIGIN, 
            radius=raio, 
            start_angle=PI/2,  # Começa no topo (0,2,0)
            angle=-PI,         # Desce até (0,-2,0) passando por (0,0,2)
            **kwargs
        )
        self.rotate(-90*DEGREES,Y_AXIS,about_point=ORIGIN)
        self.set_color(cor)
        self.set_stroke(width=espessura, opacity=opacidade)
        self.apply_depth_test()
 
 
##claude mobject
if not VMOBJECT:

    class GrandeArco(Mobject):
        """
        O que faz
        Constrói o arco de grande círculo entre dois pontos como um VMobject interno, por interpolação esférica.

        Como usar
            arco = GrandeArco(p1, p2, cor=WHITE, num_pontos=50, center=ORIGIN, espessura=LINE_SIZE)

        Parâmetros
            ponto1 (VMobject, opcional): Primeiro ponto da esfera. Padrão é None.
            ponto2 (VMobject, opcional): Segundo ponto da esfera. Padrão é None.
            cor (Color, opcional): Cor do arco. Padrão é WHITE.
            num_pontos (int, opcional): Número de amostras. Padrão é 50.
            center (np.ndarray, opcional): Centro da esfera. Padrão é ORIGIN.
            espessura (float, opcional): Largura do traço. Padrão é LINE_SIZE.
        """
        def __init__(self, ponto1, ponto2, cor=WHITE, num_pontos=50, center=None, espessura=LINE_SIZE):
            super().__init__()
            
            if center is None:
                center = ORIGIN
            else:
                center = center
                
            # Obtém as coordenadas cartesianas dos pontos na esfera
            inicio = ponto1.get_center() - center
            fim = ponto2.get_center() - center
            
            # Verifica se os pontos são iguais; se forem, não há arco a ser desenhado
            if np.allclose(inicio, fim):
                return
            
            # Calcula o ângulo entre os vetores que representam os pontos na esfera
            angulo = np.arccos(np.dot(inicio, fim) / (np.linalg.norm(inicio) * np.linalg.norm(fim)))
            
            # Lista para armazenar os pontos do arco interpolado
            pontos_arco = []
            
            # Gera pontos intermediários ao longo do arco usando interpolação esférica (Slerp)
            for t in range(num_pontos + 1):
                ponto_slerp = (
                    np.sin((1 - t / num_pontos) * angulo) * inicio +
                    np.sin((t / num_pontos) * angulo) * fim
                ) / np.sin(angulo) + center
                pontos_arco.append(ponto_slerp)
            
            # Cria uma linha poligonal conectando os pontos do arco usando VMobject
            if len(pontos_arco) > 1:
                # Cria um VMobject com os pontos
                arco_vmobject = VMobject()
                arco_vmobject.set_points_as_corners(pontos_arco)
                arco_vmobject.set_color(cor)
                arco_vmobject.set_stroke(width=espessura)
                self.add(arco_vmobject)
            self.apply_depth_test()
            
    class GrandeCirculo(Mobject):
        """
        O que faz
        Desenha um grande círculo como curva paramétrica interna e o rotaciona para alinhar com o vetor normal fornecido.

        Como usar
            gc = GrandeCirculo(vetor_normal, raio=CELESTIAL_SPHERE_RADIUS, cor=YELLOW, espessura=1.2*LINE_SIZE)

        Parâmetros
            vetor_normal (np.ndarray, opcional): Vetor normal unitário ao plano do círculo. Padrão é vetor unitário em Z.
            raio (float, opcional): Raio da esfera de referência. Padrão é CELESTIAL_SPHERE_RADIUS.
            cor (Color, opcional): Cor do traço. Padrão é YELLOW.
            espessura (float, opcional): Largura do traço. Padrão é 1.2 vezes LINE_SIZE.
        """
        def __init__(self, vetor_normal, raio=CELESTIAL_SPHERE_RADIUS, cor=YELLOW, espessura=1.2*LINE_SIZE):

            super().__init__()
            
            # Normaliza o vetor normal
            vetor_normal = vetor_normal / np.linalg.norm(vetor_normal)
            raio_corrigido = raio 
            
            # Círculo base no plano XY usando ParametricCurve
            def circulo_base(t):
                return raio_corrigido * np.array([np.cos(t), np.sin(t), 0])
            
            circulo = ParametricCurve(
                lambda t: circulo_base(t),
                t_range=[0, TAU, 0.2]
            )
            circulo.set_color(cor)
            circulo.set_stroke(width=espessura)
            
            # Calcula o ângulo de rotação entre o vetor [0, 0, 1] e o vetor normal
            z_axis = np.array([0, 0, 1])
            if not np.allclose(vetor_normal, z_axis):
                angle = np.arccos(np.dot(z_axis, vetor_normal))
                axis_of_rotation = np.cross(z_axis, vetor_normal)
                if np.linalg.norm(axis_of_rotation) > 1e-10:  # Evita divisão por zero
                    axis_of_rotation = axis_of_rotation / np.linalg.norm(axis_of_rotation)
                    circulo.rotate(angle, axis_of_rotation)
            
            self.add(circulo)
            self.apply_depth_test()

    class Paralelo(Mobject):
        """
        O que faz
        Desenha um paralelo de latitude fixa como um VMobject de linhas.

        Como usar
            paralelo = Paralelo(altitude=30, latitude=0, raio=CELESTIAL_SPHERE_RADIUS)

        Parâmetros
            altitude (float, opcional): Latitude do círculo em graus. Padrão é 0.
            latitude (float, opcional): Latitude do observador para inclinação. Padrão é 0.
            raio (float, opcional): Raio da esfera. Padrão é CELESTIAL_SPHERE_RADIUS.
            num_pontos (int, opcional): Amostragem do traço. Padrão é 400.
            espessura (float, opcional): Largura do traço. Padrão é LINE_SIZE.
            cor (Color, opcional): Cor do traço. Padrão é BLUE_C.
        """
        def __init__(self, altitude, latitude=0, raio=CELESTIAL_SPHERE_RADIUS, num_pontos=400, espessura=LINE_SIZE, cor=BLUE_C):

            super().__init__()

            # Converter latitude de graus para radianos
            latitude_rad = np.radians(altitude)

            # Gerar pontos ao longo do círculo de latitude variando a longitude (theta)
            pontos_círculo = []
            for theta in np.linspace(0, 2 * np.pi, num_pontos):
                # A posição no círculo em coordenadas cartesianas
                x = raio * np.cos(latitude_rad) * np.cos(theta)
                y = raio * np.cos(latitude_rad) * np.sin(theta)
                z = raio * np.sin(latitude_rad)
                pontos_círculo.append(np.array([x, y, z]))

            # Criar o círculo como um VMobject
            if len(pontos_círculo) > 1:
                circulo_vmobject = VMobject()
                circulo_vmobject.set_points_as_corners(pontos_círculo)
                circulo_vmobject.set_color(cor)
                circulo_vmobject.set_stroke(width=espessura)
                
                # Aplicar rotação se necessário
                if latitude != 0:
                    circulo_vmobject.rotate_about_origin((latitude - 90) * DEGREES, X_AXIS)
                
                self.add(circulo_vmobject)
            self.apply_depth_test()

    class ArcoParalelo(Mobject):
        """
        O que faz
        Constrói um arco de paralelo limitado por azimutes inicial e final.

        Como usar
            arco = ArcoParalelo(altitude=20, angulo_inicial=0, angulo_final=90, raio=CELESTIAL_SPHERE_RADIUS)

        Parâmetros
            altitude (float, opcional): Latitude do arco em graus. Padrão é 0.
            angulo_inicial (float, opcional): Azimute inicial em graus. Padrão é 0.
            angulo_final (float, opcional): Azimute final em graus. Padrão é 0.
            raio (float, opcional): Raio da esfera. Padrão é CELESTIAL_SPHERE_RADIUS.
            num_pontos (int, opcional): Pontos amostrados do arco. Padrão é 50.
            espessura (float, opcional): Largura do traço. Padrão é LINE_SIZE.
            cor (Color, opcional): Cor do arco. Padrão é WHITE.
        """
        def __init__(self, altitude, angulo_inicial, angulo_final, raio=CELESTIAL_SPHERE_RADIUS, num_pontos=50,espessura=LINE_SIZE, cor=WHITE):

            super().__init__()
            
            angulo_inicio = 90 - angulo_inicial
            angulo_fim = 90 - angulo_final
            
            # Converter graus para radianos
            latitude_rad = np.radians(altitude)
            angulo_inicio_rad = np.radians(angulo_inicio)
            angulo_fim_rad = np.radians(angulo_fim)

            # Gerar pontos apenas no intervalo especificado
            pontos_arco = []
            for theta in np.linspace(angulo_inicio_rad, angulo_fim_rad, num_pontos):
                x = raio * np.cos(latitude_rad) * np.cos(theta)
                y = raio * np.cos(latitude_rad) * np.sin(theta)
                z = raio * np.sin(latitude_rad)
                pontos_arco.append(np.array([x, y, z]))

            # Criar o arco como um VMobject
            if len(pontos_arco) > 1:
                arco_vmobject = VMobject()
                arco_vmobject.set_points_as_corners(pontos_arco)
                arco_vmobject.set_color(cor)
                self.add(arco_vmobject)
            self.apply_depth_test()
            self.set_width(espessura)

    class MeridianoLocal(Mobject):
        """
        O que faz
        Forma o meridiano local como um arco rotacionado que passa por zênite e nadir, usando Arc internamente.

        Como usar
            mer = MeridianoLocal(raio=CELESTIAL_SPHERE_RADIUS, espessura=LINE_SIZE, cor=ORANGE)

        Parâmetros
            raio (float, opcional): Raio do arco. Padrão é CELESTIAL_SPHERE_RADIUS.
            espessura (float, opcional): Largura do traço. Padrão é LINE_SIZE.
            cor (Color, opcional): Cor do traço. Padrão é ORANGE.
            opacidade (float, opcional): Opacidade do traço. Padrão é 1.
            **kwargs: Argumentos adicionais herdados de Arc.
        """
        def __init__(self, raio=CELESTIAL_SPHERE_RADIUS, espessura=LINE_SIZE, cor=ORANGE, opacidade=1, **kwargs):

            super().__init__()
            
            # Criar um arco semicircular usando Arc
            arco = Arc(
                arc_center=ORIGIN, 
                radius=raio, 
                start_angle=PI/2,  # Começa no topo
                angle=-PI,         # Desce até embaixo passando pela frente
                **kwargs
            )
            
            # Rotacionar o arco para a orientação correta
            arco.rotate(-90*DEGREES, Y_AXIS, about_point=ORIGIN)
            arco.set_color(cor)
            arco.set_stroke(width=espessura, opacity=opacidade)
            
            self.add(arco)
            self.apply_depth_test()



#DEPENDENTE DO OBSERVADOR
class EixoPolar(Line3D):
    """
    O que faz
    Cria o eixo polar inclinado de acordo com a latitude do observador, como um segmento de reta 3D.

    Como usar
        eixo = EixoPolar(latitude_graus=-23, comprimento=2*RENDER_CELESTIAL_SPHERE_RADIUS)

    Parâmetros
        latitude_graus (float, opcional): Latitude em graus para inclinação do eixo. Padrão é 0.
        comprimento (float, opcional): Comprimento total do eixo. Padrão é duas vezes RENDER_CELESTIAL_SPHERE_RADIUS.
        phi (float, opcional): Ângulo azimutal de orientação em radianos. Padrão é PI dividido por 2.
        espessura (float, opcional): Espessura do traço. Padrão é 0.03 vezes ELEMENTS_SCALE.
        cor (Color, opcional): Cor do eixo. Padrão é RED.
    """
    def __init__(self, latitude_graus, comprimento=2*RENDER_CELESTIAL_SPHERE_RADIUS, phi=PI / 2, espessura=0.03*ELEMENTS_SCALE, cor=RED):

        # Converte latitude de graus para radianos
        latitude = np.radians(latitude_graus)
        
        # Define os pontos superior e inferior do eixo no espaço 3D
        ponto_superior = np.array([ 
            comprimento / 2 * np.cos(phi) * np.cos(latitude),  # X
            comprimento / 2 * np.sin(phi) * np.cos(latitude),  # Y
            comprimento / 2 * np.sin(latitude)                # Z
        ])
        ponto_inferior = np.array([
            -comprimento / 2 * np.cos(phi) * np.cos(latitude),  # X
            -comprimento / 2 * np.sin(phi) * np.cos(latitude),  # Y
            -comprimento / 2 * np.sin(latitude)                # Z
        ])
        super().__init__(ponto_superior,ponto_inferior,width=espessura)
        self.set_color(cor)

class Equador(GrandeCirculo):
    """
    O que faz
    Desenha o círculo máximo correspondente ao equador celeste para uma latitude dada do observador.

    Como usar
        eq = Equador(latitude=-23, cor=YELLOW_D, raio=RENDER_CELESTIAL_SPHERE_RADIUS)

    Parâmetros
        latitude (float, opcional): Latitude do observador em graus. Padrão é 0.
        cor (Color, opcional): Cor do traço. Padrão é YELLOW_D.
        raio (float, opcional): Raio da esfera. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
    """
    def __init__(self, latitude, cor=YELLOW_D, raio=RENDER_CELESTIAL_SPHERE_RADIUS, **kwargs):

        # Converte a latitude para radianos
        latitude_rad = np.radians(latitude)
        
        
        
        # Define o vetor normal ao plano do equador
        vetor_normal = np.array([0, np.cos(latitude_rad), np.sin(latitude_rad)])
        
        # Chama o construtor da classe pai para criar o grande círculo
        super().__init__(vetor_normal=vetor_normal, raio=raio, cor=cor, **kwargs)

class Grade(VGroup):
    """
    O que faz
    Desenha uma grade esférica com paralelos de declinação e meridianos de ascensão reta para referência visual.

    Como usar
        grade = Grade(latitude=90, raio=RENDER_CELESTIAL_SPHERE_RADIUS, n_ar=24)

    Parâmetros
        latitude (float, opcional): Latitude para inclinar a grade. Padrão é 90.
        raio (float, opcional): Raio da esfera de referência. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        cor_dec (Color, opcional): Cor dos paralelos de declinação. Padrão é GREY_A.
        cor_ar (Color, opcional): Cor dos meridianos de ascensão reta. Padrão é GREY_A.
        opacidade (float, opcional): Opacidade das linhas. Padrão é 0.8.
        espessura (float, opcional): Largura das linhas. Padrão é 0.5 vezes LINE_SIZE.
        n_ar (int, opcional): Número de divisões de ascensão reta. Padrão é 24.
        **kwargs: Argumentos adicionais herdados de VGroup.
    """
    def __init__(self, latitude=90, raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor_dec=GREY_A,cor_ar=GREY_A, opacidade=0.8,espessura=0.5*LINE_SIZE, n_ar=24,**kwargs):
        super().__init__(**kwargs)
        self.raio = raio
        paralelos = VGroup()
        # Declinações de -80° a +80° de 10 em 10 (17 no total)
        for dec_graus in range(-80, 81, 10):
            ang = np.radians(dec_graus)  # converte para radianos
            r = raio * np.cos(ang)
            z = raio * np.sin(ang)
            paralelo = Circle(radius=r).set_color(cor_dec).set_stroke(width=espessura,opacity=opacidade)
            paralelo.shift(OUT * z)
            paralelos.add(paralelo)
        meridianos = VGroup()
        # 24 linhas de ascensão reta = 12 círculos de meridiano
        for j in range(n_ar):
            theta = TAU * j / n_ar
            arco = Arc(
                radius=raio,
                start_angle=-PI/2,
                angle=PI,
                stroke_opacity=opacidade
            ).set_color(cor_ar).set_stroke(width=espessura,opacity=opacidade)
            arco.rotate(PI/2, axis=RIGHT, about_point=ORIGIN)
            arco.rotate(theta, axis=OUT, about_point=ORIGIN)
            meridianos.add(arco)
        self.add(paralelos,meridianos)
        if latitude != 90:
            self.rotate(-90+latitude,X_AXIS,about_point=ORIGIN)
        self.apply_depth_test()

class GradeMesh(SurfaceMesh):
    """
    O que faz
    Cria uma malha de grade sobre a esfera para visualização de linhas estruturais com controle de resolução.

    Como usar
        mesh = GradeMesh(latitude=90, raio=RENDER_CELESTIAL_SPHERE_RADIUS, n_dec=19, n_ra=25)

    Parâmetros
        latitude (float, opcional): Latitude para inclinar a malha. Padrão é 90.
        raio (float, opcional): Raio da esfera de referência. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS.
        cor (Color, opcional): Cor das linhas da malha. Padrão é BLUE_C.
        n_dec (int, opcional): Resolução em declinação. Padrão é 19.
        n_ra (int, opcional): Resolução em ascensão reta. Padrão é 25.
        opacidade (float, opcional): Opacidade das linhas. Padrão é 0.5.
        espessura (float, opcional): Espessura das linhas. Padrão é 0.5 vezes LINE_SIZE.
        **kwargs: Argumentos adicionais herdados de SurfaceMesh.

    Observações
        Usa uma esfera base para gerar a malha e aplica rotação se a latitude for diferente de noventa.
    """
    def __init__(self, latitude=90,raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=BLUE_C, n_dec=19,n_ra=25, opacidade=0.5, espessura=0.5*LINE_SIZE,**kwargs):
        self.raio = raio

        def uv_func(u, v):
            return np.array([
                raio * np.cos(u) * v,
                raio * np.sin(u) * v,
                0
            ])

        super().__init__(
            Sphere(radius=raio),
            stroke_width=espessura,
            stroke_color=cor,
            resolution=(n_ra,n_dec),
            # u_range=[0, TAU],
            # v_range=[0, 1],
            **kwargs
        )
        self.set_fill(opacity=0)
        if latitude != 90:
            self.rotate(-90+latitude,X_AXIS,about_point=ORIGIN)
        self.apply_depth_test()



#ELEMENTOS EXTRA

class CirculoTangente(Polygon):
    """
    O que faz
    Constrói um círculo no plano tangente à esfera no ponto fornecido, usando o ponto como centro local.

    Como usar
        ct = CirculoTangente(ponto, raio_circulo=0.2, cor=WHITE, espessura=3)

    Parâmetros
        ponto (VMobject, opcional): Ponto na superfície da esfera que define o plano tangente. Padrão é None.
        raio_circulo (float, opcional): Raio do círculo no plano tangente. Padrão é 0.2.
        cor (Color, opcional): Cor do contorno. Padrão é WHITE.
        espessura (float, opcional): Largura do contorno. Padrão é 3.

    Observações
        Gera dois vetores tangentes ortogonais para parametrizar a circunferência no plano tangente.
    """
    def __init__(self, ponto: VMobject, raio_circulo=0.2, cor=WHITE, espessura=3):

        # Obter a posição do ponto p
        p_pos = ponto.get_center()

        # O vetor radial do ponto (do centro da esfera ao ponto 'p')
        vetor_radial = p_pos / np.linalg.norm(p_pos) 

        # O vetor normal ao plano tangente (que é o vetor radial)
        vetor_normal = vetor_radial

        # Gerando um sistema de coordenadas ortogonal (vetores tangentes ao plano)
        # Vamos usar dois vetores ortogonais ao vetor normal (usando o produto vetorial)
        # Geramos um vetor qualquer para ser ortogonal ao normal e depois normalizamos.
        vetor_qualquer = np.array([1, 0, 0]) if vetor_normal[0] != 1 else np.array([0, 1, 0])
        vetor_tangente1 = np.cross(vetor_normal, vetor_qualquer)
        vetor_tangente1 /= np.linalg.norm(vetor_tangente1)  # Normaliza

        # O segundo vetor tangente é ortogonal ao primeiro e ao normal (produto vetorial)
        vetor_tangente2 = np.cross(vetor_normal, vetor_tangente1)
        vetor_tangente2 /= np.linalg.norm(vetor_tangente2)  # Normaliza

        # Gerando o círculo no plano tangente usando uma parametrização circular
        pontos_circulo = []
        for angulo in np.linspace(0, 2 * np.pi, 100):
            ponto_no_circulo = p_pos + raio_circulo * (np.cos(angulo) * vetor_tangente1 + np.sin(angulo) * vetor_tangente2)
            pontos_circulo.append(ponto_no_circulo)
        
        # Criando o círculo com o Polygon
        super().__init__(*pontos_circulo, color=cor, stroke_width=espessura)
        self.apply_depth_test()

class AnguloEsferico(VGroup):
    """
    O que faz
    Desenha um arco indicador do ângulo esférico interno em um ponto da esfera entre dois grandes arcos, com possibilidade de rótulo matemático.

    Como usar
        ang = AnguloEsferico(p, p1, p2, raio_circulo=0.3, cor=WHITE, math_label="\\alpha")

    Parâmetros
        p (VMobject, opcional): Ponto do vértice do ângulo na esfera. Padrão é None.
        p1 (VMobject, opcional): Primeiro ponto para definição do plano. Padrão é None.
        p2 (VMobject, opcional): Segundo ponto para definição do plano. Padrão é None.
        raio_circulo (float, opcional): Raio do arco no plano tangente. Padrão é 0.3 vezes ELEMENTS_SCALE.
        cor (Color, opcional): Cor do arco. Padrão é WHITE.
        math_label (str, opcional): Texto matemático LaTeX para exibir ao lado do arco. Padrão é None.
        reducing_factor (float, opcional): Fator de recuo em relação à normal, para evitar z fighting. Padrão é 50.
        cor_do_texto (Color, opcional): Cor do rótulo. Padrão é WHITE.
        label_distance (float, opcional): Distância do rótulo ao arco. Padrão é 1.
        tamanho_da_fonte (int, opcional): Tamanho da fonte em pontos. Padrão é 30.
        num_pontos (int, opcional): Número de pontos do arco. Padrão é 30.
        center (np.ndarray, opcional): Centro geométrico da esfera caso deslocada. Padrão é None.
        espessura (float, opcional): Espessura do arco. Padrão é 2.

    Retorno
        VGroup: Grupo contendo o arco e, se definido, o rótulo.

    Observações
        O valor do ângulo é calculado com a função angulo_diedro_de_vetores_cartesianos e pode ser obtido por valor_angulo.
    """
    def __init__(self, p: VMobject, p1: VMobject, p2: VMobject, raio_circulo=0.3*ELEMENTS_SCALE, cor=WHITE, math_label = None,reducing_factor = 50, cor_do_texto=WHITE, label_distance=1,tamanho_da_fonte=30, num_pontos=30, center=None,espessura=2):

        
        if center is None:
            center = ORIGIN
        else:
            center = center
        self.raio_circulo = raio_circulo
        self.cor = cor
        self.num_pontos = num_pontos
        self.espesura = espessura
        # Chama o mÃ©todo da classe base (CirculoTangente) para definir a posiÃ§Ã£o e os vetores tangentes
        p_pos = p.get_center() - center
        self.p_pos = p_pos
        
        self.p1 = p1
        self.p2 = p2
        p1_pos = self.p1.get_center() - center
        p2_pos = self.p2.get_center() - center
        normal_vector = p_pos  # Vetor normal ao plano tangente

        # Criar dois vetores tangentes ao plano
        vetor_qualquer = np.array(p1_pos/np.linalg.norm(p1_pos))
        vetor_tangente1 = np.cross(normal_vector, vetor_qualquer)
        vetor_tangente1 /= np.linalg.norm(vetor_tangente1)
        
        vetor_tangente2 = np.cross(normal_vector, vetor_tangente1)
        vetor_tangente2 /= np.linalg.norm(vetor_tangente2)

        # Converter Ã¢ngulos para radianos
        self.angulo_final = np.radians(angulo_diedro_de_vetores_cartesianos(p1_pos, p2_pos, p_pos))
        
        self._valor_angulo = float(np.round(np.degrees(self.angulo_final),2))
        # Gerar pontos do arco
        pontos_arco = [
            center + 1.01*p_pos + raio_circulo * (np.cos(angulo) * -vetor_tangente2 + np.sin(angulo) * -vetor_tangente1) - normal_vector / reducing_factor
            for angulo in np.linspace(np.radians(0), self.angulo_final, num_pontos)
        ]
        super().__init__()
        # Criar o arco como um VMobject (para suavizaÃ§Ã£o)
        arco2 = VMobject().set_stroke(width=espessura,color=cor)  # Chama a inicializaÃ§Ã£o de CirculoTangente
        arco2.set_points_as_corners(pontos_arco)
        if math_label != None:
            texto = Tex(fr"{math_label}", font_size = tamanho_da_fonte).move_to(p_pos + label_distance*1.7*raio_circulo * (np.cos(self.angulo_final/2) * -vetor_tangente2 + np.sin(self.angulo_final/2) * -vetor_tangente1)).set_color(cor_do_texto)
            self.add(texto)
        self.add(arco2)
        self.apply_depth_test()
    
    def valor_angulo(self) -> float:  # Agora Ã© acessado como um atributo
        """
        O que faz
        Retorna o valor atual do ângulo esférico interno no vértice definido.

        Como usar
            valor = angulo.valor_angulo()

        Parâmetros
            Não se aplica.

        Retorno
            float: Ângulo em graus calculado entre os planos definidos.
        """
        return angulo_diedro_de_vetores_cartesianos(self.p1.get_center(), self.p2.get_center(), self.p_pos)
   
    def atualizar_angulo(self, p: VMobject, p1: VMobject, p2: VMobject):
        """
        O que faz
        Recalcula e atualiza o arco do ângulo esférico após mover o vértice ou os pontos de referência.

        Como usar
            angulo.atualizar_angulo(novo_p, novo_p1, novo_p2)

        Parâmetros
            p (VMobject, opcional): Novo ponto do vértice. Padrão é None.
            p1 (VMobject, opcional): Novo primeiro ponto de referência. Padrão é None.
            p2 (VMobject, opcional): Novo segundo ponto de referência. Padrão é None.

        Retorno
            None: Atualiza os pontos internos do objeto.
        """
        p_pos = p.get_center()
        self.p_pos = p_pos
        
        self.p1 = p1
        self.p2 = p2
        
        normal_vector = p_pos / np.linalg.norm(p_pos)  # Vetor normal ao plano tangente

        # Criar dois vetores tangentes ao plano
        vetor_qualquer = np.array(p1.get_center()/np.linalg.norm(p1.get_center()))
        vetor_tangente1 = np.cross(normal_vector, vetor_qualquer)
        vetor_tangente1 /= np.linalg.norm(vetor_tangente1)
        
        vetor_tangente2 = np.cross(normal_vector, vetor_tangente1)
        vetor_tangente2 /= np.linalg.norm(vetor_tangente2)

        # Converter Ã¢ngulos para radianos
        self.angulo_final = np.radians(angulo_diedro_de_vetores_cartesianos(p1.get_center(), p2.get_center(), p_pos))
        
        self._valor_angulo = float(np.round(np.degrees(self.angulo_final),2))
        # Gerar pontos do arco
        pontos_arco = [
            p_pos + self.raio_circulo * (np.cos(angulo) * -vetor_tangente2 + np.sin(angulo) * -vetor_tangente1) - normal_vector / 30
            for angulo in np.linspace(np.radians(0), self.angulo_final, self.num_pontos)
        ]
        
        # Criar o arco como um VMobject (para suavizaÃ§Ã£o)
         # Chama a inicializaÃ§Ã£o de CirculoTangente
        self.set_points_as_corners(pontos_arco)
     
class RegiaoLatitudinal(Surface):
    """
    O que faz
    Cria uma faixa esférica delimitada por duas latitudes, útil para destacar zonas como trópicos ou círculos polares.

    Como usar
        faixa = RegiaoLatitudinal(lat_ini=-23.5, lat_fim=23.5, raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=YELLOW, opacidade=0.5)

    Parâmetros
        lat_ini (float, opcional): Latitude inicial em graus. Padrão é 0.
        lat_fim (float, opcional): Latitude final em graus. Padrão é 0.
        raio (float, opcional): Raio da esfera. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS mais 0.02.
        cor (Color, opcional): Cor de preenchimento. Padrão é YELLOW.
        opacidade (float, opcional): Opacidade do preenchimento. Padrão é 0.5.
        resolucao (tuple, opcional): Resolução da malha u e v. Padrão é 40, 40.
        **kwargs: Argumentos adicionais herdados de Surface.

    Observações
        A parametrização percorre a longitude completa e interpola a latitude entre os limites informados.
    """ 
    def __init__(self, lat_ini, lat_fim, raio=RENDER_CELESTIAL_SPHERE_RADIUS+0.02, cor=YELLOW, opacidade=0.5, resolucao=(40, 40), **kwargs):

        lat_ini = np.radians(lat_ini)
        lat_fim = np.radians(lat_fim)

        self.raio= raio
        self.lat_ini,self.lat_fim=lat_ini, lat_fim
      

        super().__init__(u_range=(0, 1), v_range=(0, 1), resolution=resolucao, **kwargs)
        self.set_color(cor)  # Define a cor da região
        self.set_opacity(opacidade) 
        
    def uv_func(self, u, v):
        """ Converte coordenadas paramétricas (u, v) em pontos 3D na esfera """
        theta = interpolate(self.lat_ini, self.lat_fim, u)  # Latitude interpolada
        phi = interpolate(0, 2 * np.pi, v)        # Longitude completa

        x = self.raio * np.cos(theta) * np.cos(phi)
        y = self.raio * np.cos(theta) * np.sin(phi)
        z = self.raio * np.sin(theta)

        return np.array([x, y, z])
        
class RegiaoMeridional(Surface):
    """
    O que faz
    Cria uma faixa esférica delimitada por duas longitudes, útil para destacar setores específicos como intervalos de azimute.

    Como usar
        setor = RegiaoMeridional(lon_ini=0, lon_fim=60, raio=RENDER_CELESTIAL_SPHERE_RADIUS, cor=GREEN, opacidade=0.5)

    Parâmetros
        lon_ini (float, opcional): Longitude inicial em graus. Padrão é 0.
        lon_fim (float, opcional): Longitude final em graus. Padrão é 0.
        raio (float, opcional): Raio da esfera. Padrão é RENDER_CELESTIAL_SPHERE_RADIUS mais 0.02.
        cor (Color, opcional): Cor de preenchimento. Padrão é GREEN.
        opacidade (float, opcional): Opacidade do preenchimento. Padrão é 0.5.
        resolucao (tuple, opcional): Resolução da malha u e v. Padrão é 40, 40.
        **kwargs: Argumentos adicionais herdados de Surface.

    Observações
        A latitude é varrida de menos noventa a mais noventa, enquanto a longitude é interpolada entre os limites definidos.
    """
    def __init__(self, lon_ini, lon_fim, raio=RENDER_CELESTIAL_SPHERE_RADIUS+0.02, cor=GREEN, opacidade=0.5, resolucao=(40, 40), **kwargs):
        
        lon_ini = np.radians(lon_ini)
        lon_fim = np.radians(lon_fim)
        self.raio = raio
        self.lon_ini = lon_ini
        self.lon_fim = lon_fim

        

        super().__init__(u_range=(0, 1), v_range=(0, 1), resolution=resolucao, **kwargs)
        
        self.set_color(cor)  # Define a cor da região
        self.set_opacity(opacidade)
    def uv_func(self, u, v):
            """ Converte coordenadas paramétricas (u, v) em pontos 3D na esfera """
            phi = interpolate(self.lon_ini, self.lon_fim, u)  # Longitude interpolada
            theta = interpolate(-np.pi / 2, np.pi / 2, v)  # Latitude total

            x = self.raio * np.cos(theta) * np.cos(phi)
            y = self.raio * np.cos(theta) * np.sin(phi)
            z = self.raio * np.sin(theta)

            return np.array([x, y, z])

class LinhaVMobject(VMobject):
    """
    O que faz
    Cria uma linha entre dois pontos arbitrários no espaço como um VMobject simples.

    Como usar
        linha = LinhaVMobject(ponto_inicial, ponto_final, cor=WHITE, espessura=2)

    Parâmetros
        ponto_inicial (np.ndarray, opcional): Coordenadas do ponto inicial. Padrão é None.
        ponto_final (np.ndarray, opcional): Coordenadas do ponto final. Padrão é None.
        cor (Color, opcional): Cor do traço. Padrão é WHITE.
        espessura (float, opcional): Largura do traço. Padrão é 2.
        **kwargs: Argumentos adicionais herdados de VMobject.
    """
    def __init__(self, ponto_inicial, ponto_final, cor=WHITE, espessura=2, **kwargs):
        super().__init__(**kwargs)
        
        self.set_points_as_corners([ponto_inicial, ponto_final])
        self.set_color(cor)
        self.set_stroke(width=espessura)

class Seta(Arrow):
    """
    O que faz
    Cria uma seta tangente à esfera no ponto dado, orientada por um ângulo relativo ao norte local.

    Como usar
        seta = Seta(ponto, angulo=45, comprimento=0.5, cor=WHITE, espessura=2)

    Parâmetros
        ponto (VMobject, opcional): Ponto de origem da seta sobre a esfera. Padrão é None.
        angulo (float, opcional): Ângulo em graus em relação à direção norte local. Padrão é 0.
        comprimento (float, opcional): Comprimento da seta. Padrão é 0.5.
        cor (Color, opcional): Cor da seta. Padrão é WHITE.
        espessura (float, opcional): Largura do traço. Padrão é 2.

    Observações
        Usa a base tangente local formada por dois vetores ortogonais ao vetor radial para orientar a seta.
    """
    def __init__(self, ponto, angulo: float, comprimento=0.5, cor=WHITE, espessura=2):

        # Obter a posição do ponto
        p_pos = ponto.get_center()

        # O vetor radial do ponto (do centro da esfera ao ponto 'p')
        vetor_radial = p_pos / np.linalg.norm(p_pos)

        # Gerando um sistema de coordenadas ortogonal
        vetor_qualquer = np.array([0, 0, 1])
        vetor_tangente1 = np.cross(vetor_radial, vetor_qualquer)
        vetor_tangente1 /= np.linalg.norm(vetor_tangente1)
        
        vetor_tangente2 = np.cross(vetor_radial, vetor_tangente1)
        vetor_tangente2 /= np.linalg.norm(vetor_tangente2)
        
        # Criar o vetor da seta inclinada do ângulo dado
        direcao_seta = np.cos(np.radians(angulo)) * -vetor_tangente2 + np.sin(np.radians(angulo)) * vetor_tangente1
        direcao_seta *= comprimento  # Ajusta o comprimento

        # Ponto final da seta
        ponto_final = p_pos + direcao_seta
        
        # Criando a seta com Arrow
        super().__init__(start=p_pos, end=ponto_final, color=cor, stroke_width=espessura, buff=0)

class SetaGrandeCirculo(Seta):
    """
    O que faz
    Cria uma seta tangente ao grande círculo especificado no ponto dado, apontando ao longo da direção local do círculo.

    Como usar
        seta_gc = SetaGrandeCirculo(ponto, grande_circulo, comprimento=0.5, cor=WHITE, espessura=2, inverter=False)

    Parâmetros
        ponto (VMobject, opcional): Ponto na esfera onde a seta nasce. Padrão é None.
        grande_circulo (VGroup, opcional): Objeto contendo o grande círculo desenhado. Padrão é None.
        comprimento (float, opcional): Comprimento da seta. Padrão é 0.5.
        cor (Color, opcional): Cor da seta. Padrão é WHITE.
        espessura (float, opcional): Largura do traço. Padrão é 2.
        inverter (bool, opcional): Indica se a direção deve ser invertida. Padrão é False.
        granularidade (int, opcional): Redução de pontos para busca de tangente. Padrão é 100.

    Observações
        Seleciona o ponto mais próximo no grande círculo e usa o vizinho para estimar a direção tangente.
    """
    def __init__(self, ponto, grande_circulo: VGroup, comprimento=0.5, cor=WHITE, espessura=2, inverter=False, granularidade=100):

        # Obter a posição do ponto na esfera
        p_pos = ponto.get_center()

        # Coletar os pontos de todos os componentes do VGroup
        pontos_circulo = np.vstack([obj.points for obj in grande_circulo if hasattr(obj, "points") and len(obj.points) > 0])

        # Verificar se conseguimos pontos válidos
        if len(pontos_circulo) == 0:
            raise ValueError("O grande círculo (VGroup) não contém pontos!")

        # Reduzindo a quantidade de pontos para melhorar o desempenho
        pontos_circulo = pontos_circulo[:: max(1, len(pontos_circulo) // granularidade)]

        # Encontrar o ponto mais próximo no grande círculo
        ponto_mais_proximo = min(pontos_circulo, key=lambda pt: np.linalg.norm(pt - p_pos))

        # Encontrar o índice do ponto mais próximo
        idx_ponto = np.where((pontos_circulo == ponto_mais_proximo).all(axis=1))[0]
        if len(idx_ponto) == 0:
            raise ValueError("Não foi possível encontrar o índice do ponto mais próximo!")
        idx_ponto = idx_ponto[0]

        # Escolher um ponto vizinho para calcular a direção tangente
        ponto_vizinho = pontos_circulo[(idx_ponto + 1) % len(pontos_circulo)]

        # Vetor tangente ao grande círculo no ponto dado
        vetor_tangente = ponto_vizinho - ponto_mais_proximo

        # Se o vetor for zero, lançar um erro
        if np.linalg.norm(vetor_tangente) < 1e-6:
            raise ValueError("Vetor tangente muito pequeno ou inválido!")

        # Normaliza o vetor
        vetor_tangente /= np.linalg.norm(vetor_tangente)

        # Inverter a direção da seta se necessário
        if inverter:
            vetor_tangente *= -1

        # Criando a seta apontando na direção do grande círculo
        super().__init__(ponto, angulo=0, comprimento=comprimento, cor=cor, espessura=espessura)

        # Atualizando a direção da seta com o vetor tangente calculado
        self.put_start_and_end_on(p_pos, p_pos + vetor_tangente * comprimento)

        # Ajustando para melhor renderização 3D
        self.set_shade_in_3d(True)

class PontosCardeais(VGroup):
    """
    Cria um objeto que representa os pontos cardeais (Norte, Sul, Leste, Oeste) em torno de um círculo.

    Parâmetros:
        raio (float): Raio do círculo (padrão: 2).
        color (Color): Cor dos pontos cardeais e do círculo (padrão: WHITE).
        font_size (int): Tamanho da fonte dos pontos cardeais (padrão: 40).
        weight (str): Peso da fonte (padrão: BOLD).
    """
    def __init__(self, raio=1*RENDER_CELESTIAL_SPHERE_RADIUS, color=WHITE, font_size=10, weight=BOLD):
        super().__init__()

        # Cria o círculo central para os pontos cardeais
        mobject = Circle(radius=raio, color=color)

        # Cria os textos para os pontos cardeais (Norte, Sul, Leste, Oeste)
        north = Text("N", font_size=font_size, weight=weight, color=color).next_to(mobject.get_top(), UP/4,buff=0).set_z_index(2)
        south = Text('S', font_size=font_size, weight=weight, color=color).next_to(mobject.get_bottom(), DOWN/4,buff=0).set_z_index(2)
        west = Text("O", font_size=font_size, weight=weight, color=color).next_to(mobject.get_left(), LEFT/4,buff=0).set_z_index(2)
        east = Text("L", font_size=font_size, weight=weight, color=color).next_to(mobject.get_right(), RIGHT/4,buff=0).set_z_index(2)

        # Adiciona os pontos cardeais ao grupo
        self.add(north, south, west, east)

#MARCADORES

class MarcadorAngulo(VGroup):
    """
    O que faz
    Marca um ângulo no espaço entre duas linhas radiais que partem de uma origem até dois pontos, com arco opcional e rótulos.

    Como usar
        marcador = MarcadorAngulo(p1, p2, raio_arco=0.5, barra=True, cor_linha=WHITE, cor_arco=WHITE)

    Parâmetros
        ponto1 (VMobject ou np.ndarray, opcional): Primeiro ponto do ângulo ou vetor posição. Padrão é None.
        ponto2 (VMobject ou np.ndarray, opcional): Segundo ponto do ângulo ou vetor posição. Padrão é None.
        raio_arco (float, opcional): Raio do arco que representa o ângulo. Padrão é 0.5.
        Origem (VMobject ou np.ndarray, opcional): Objeto ou vetor da origem das linhas. Padrão é None para ORIGIN.
        barra (bool, opcional): Indica se desenha as barras radiais. Padrão é True.
        cor_linha (Color, opcional): Cor das barras radiais. Padrão é WHITE.
        espessura_linha (float, opcional): Espessura das barras radiais. Padrão é 2.
        cor_arco (Color, opcional): Cor do arco do ângulo. Padrão é WHITE.
        cor_texto (Color, opcional): Cor do texto do rótulo. Padrão é WHITE.
        espessura_arco (float, opcional): Espessura do arco. Padrão é 3.
        label_2d (bool, opcional): Indica rótulo 2D. Padrão é False.
        label_3d (bool, opcional): Indica rótulo 3D. Padrão é False.
        math_label_2d (str, opcional): Texto LaTeX do rótulo 2D. Padrão é None.
        math_label_3d (str, opcional): Texto LaTeX do rótulo 3D. Padrão é None.
        label_distance (float, opcional): Distância do rótulo ao arco. Padrão é 1.
        tamanho_da_fonte (int, opcional): Tamanho do rótulo. Padrão é 20.
        arco (bool, opcional): Indica se desenha o arco. Padrão é True.
        **kwargs: Argumentos adicionais herdados de VGroup.
    """
    def __init__(self, ponto1:VMobject, ponto2:VMobject,raio_arco=0.5, Origem = None , barra=True, cor_linha=WHITE, espessura_linha=2, cor_arco=WHITE, cor_texto=WHITE, espessura_arco=3, label_2d=False,label_3d=False, math_label_2d=None, math_label_3d=None, label_distance=1, tamanho_da_fonte=20, arco=True, **kwargs):

        super().__init__(**kwargs)
        if Origem == None:
            Origem = ORIGIN
        else: 
            Origem = Origem.get_center()
        # Obtém as posições dos dois pontos no espaço
       
        
        if isinstance(ponto1, VMobject)or isinstance(ponto1, P):
            ponto1_pos = ponto1.get_center()
        else:
            ponto1_pos = ponto1
        if isinstance(ponto2, VMobject)or isinstance(ponto2, P):
            ponto2_pos = ponto2.get_center()
        else:
            ponto2_pos = ponto2
        
        # Cria as duas linhas radiais (do centro aos pontos)
        segmentos_linha1 = LinhaVMobject(Origem,ponto1_pos,cor_linha,espessura_linha)
        segmentos_linha2 = LinhaVMobject(Origem,ponto2_pos,cor_linha,espessura_linha)

    
        # Obtém os vetores das duas linhas
        linha1_vector = ponto1_pos - Origem
        linha2_vector = ponto2_pos - Origem

        # Cálculo do ângulo entre os dois vetores usando o produto escalar
        cos_theta = np.dot(linha1_vector, linha2_vector) / (np.linalg.norm(linha1_vector) * np.linalg.norm(linha2_vector))
        angulo = np.arccos(cos_theta)  # Em radianos

        # Agora, criamos um arco entre os dois vetores no plano formado pelo produto vetorial
        # Vetor perpendicular ao plano formado pelos vetores linha1_vector e linha2_vector
        produto_vetorial = np.cross(linha1_vector, linha2_vector)
        produto_vetorial /= np.linalg.norm(produto_vetorial) 
        
        produto_vetorial2 = np.cross(linha1_vector, produto_vetorial)
        produto_vetorial2 /= np.linalg.norm(produto_vetorial2)

        num_pontos = 50
        pontos_arco = [
                Origem + np.cos(angulo) * linha1_vector / np.linalg.norm(linha1_vector) * 0.5 *raio_arco + np.sin(angulo) * -produto_vetorial2 / np.linalg.norm(produto_vetorial2) * 0.5*raio_arco
                for angulo in np.linspace(np.radians(0), angulo, num_pontos)
            ]
            
        arco2 = VMobject().set_stroke(width=espessura_arco,color=cor_arco)  # Chama a inicializaÃ§Ã£o de CirculoTangente
        arco2.set_points_as_corners(pontos_arco)
        
        if label_2d:
            texto = Tex(fr"{np.round(np.degrees(angulo),2)}^\circ", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
            texto.rotate(PI/2, X_AXIS).rotate(PI/2, Z_AXIS)
            self.add(texto)
        if label_3d:
            texto = Tex(fr"{np.round(np.degrees(angulo),2)}^\circ", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
            self.add(texto)
        if math_label_2d != None:
            texto = Tex(fr"{math_label_2d}", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center()).set_color(cor_texto)
            texto.rotate(PI/2, X_AXIS).rotate(PI/2, Z_AXIS)
            self.add(texto)
        if math_label_3d != None:
            texto = Tex(fr"{math_label_3d}", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center()).set_color(cor_texto)
            self.add(texto)
            
      

        # Adiciona as linhas e o arco ao VGroup
        if barra:
            self.add(segmentos_linha1, segmentos_linha2)
        if arco:
            self.add(arco2)

class MarcadorAltura(VGroup):
    """
    O que faz
    Marca a altura de um ponto acima do plano XY por um arco no plano vertical definido pela projeção do ponto no plano.

    Como usar
        marcador = MarcadorAltura(ponto, cor_linha=WHITE, cor_arco=WHITE)

    Parâmetros
        ponto (VMobject, opcional): Ponto cuja altura será marcada. Padrão é None.
        Origem (VMobject ou np.ndarray, opcional): Origem das linhas. Padrão é ORIGIN.
        barra (bool, opcional): Indica se desenha a barra radial. Padrão é True.
        cor_linha (Color, opcional): Cor da barra radial. Padrão é WHITE.
        espessura_linha (float, opcional): Espessura da barra radial. Padrão é 2.
        cor_arco (Color, opcional): Cor do arco do ângulo. Padrão é WHITE.
        cor_texto (Color, opcional): Cor dos rótulos. Padrão é WHITE.
        espessura_arco (float, opcional): Espessura do arco. Padrão é 3.
        label_2d (bool, opcional): Indica rótulo 2D. Padrão é False.
        label_3d (bool, opcional): Indica rótulo 3D. Padrão é False.
        math_label_2d (str, opcional): Texto do rótulo 2D. Padrão é None.
        math_label_3d (str, opcional): Texto do rótulo 3D. Padrão é None.
        label_distance (float, opcional): Distância do rótulo ao arco. Padrão é 1.
        tamanho_da_fonte (int, opcional): Tamanho da fonte. Padrão é 20.
        arco (bool, opcional): Indica se desenha o arco de altura. Padrão é True.
        **kwargs: Argumentos adicionais herdados de VGroup.
    """
    def __init__(self, ponto:VMobject, Origem = None , barra=True, cor_linha=WHITE, espessura_linha=2, cor_arco=WHITE, cor_texto=WHITE, espessura_arco=3, label_2d=False,label_3d=False, math_label_2d=None, math_label_3d=None, label_distance=1, tamanho_da_fonte=20, arco=True, **kwargs):

        super().__init__(**kwargs)
        if Origem == None:
            Origem = ORIGIN
        else: 
            Origem = Origem.get_center()
        # Obtém as posições dos dois pontos no espaço
        ponto1_pos = ponto.get_center()
        ponto2_pos = [ponto1_pos[0],ponto1_pos[1],0]
        class LinhaVMobject(VMobject):
            def __init__(self, ponto_inicial, ponto_final, cor=WHITE, espessura=2, **kwargs):
                super().__init__(**kwargs)
                self.set_points_as_corners([ponto_inicial, ponto_final])
                self.set_color(cor)
                self.set_stroke(width=espessura)
        # Cria as duas linhas radiais (do centro aos pontos)
        segmentos_linha1 = LinhaVMobject(Origem,ponto1_pos,cor_linha,espessura_linha)

    
        # Obtém os vetores das duas linhas
        linha1_vector = ponto1_pos - Origem
        linha2_vector = ponto2_pos - Origem

        # Cálculo do ângulo entre os dois vetores usando o produto escalar
        cos_theta = np.dot(linha1_vector, linha2_vector) / (np.linalg.norm(linha1_vector) * np.linalg.norm(linha2_vector))
        angulo = np.arccos(cos_theta)  # Em radianos

        # Agora, criamos um arco entre os dois vetores no plano formado pelo produto vetorial
        # Vetor perpendicular ao plano formado pelos vetores linha1_vector e linha2_vector
        produto_vetorial = np.cross(linha1_vector, linha2_vector)
        produto_vetorial /= np.linalg.norm(produto_vetorial) 
        
        produto_vetorial2 = np.cross(linha1_vector, produto_vetorial)
        produto_vetorial2 /= np.linalg.norm(produto_vetorial2)

        num_pontos = 50
        pontos_arco = [
                Origem + np.cos(angulo) * linha1_vector / np.linalg.norm(linha1_vector) * 0.5 + np.sin(angulo) * -produto_vetorial2 / np.linalg.norm(produto_vetorial2) * 0.5
                for angulo in np.linspace(np.radians(0), angulo, num_pontos)
            ]
            
        arco2 = VMobject().set_stroke(width=espessura_arco,color=cor_arco)  # Chama a inicializaÃ§Ã£o de CirculoTangente
        arco2.set_points_as_corners(pontos_arco)
        
        if label_2d:
            texto = Tex(fr"{np.round(np.degrees(angulo),2)}^\circ", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
            texto.rotate(PI/2, X_AXIS).rotate(PI/2, Z_AXIS)
            self.add(texto)
        if label_3d:
            texto = Tex(fr"{np.round(np.degrees(angulo),2)}^\circ", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
            self.add(texto)
        if math_label_2d != None:
            texto = Tex(fr"{math_label_2d}", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
            texto.rotate(PI/2, X_AXIS).rotate(PI/2, Z_AXIS)
            self.add(texto)
        if math_label_3d != None:
            texto = Tex(fr"{math_label_3d}", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
            self.add(texto)
            
      

        # Adiciona as linhas e o arco ao VGroup
        if barra:
            self.add(segmentos_linha1)
        if arco:
            self.add(arco2)

class MarcadorLatitude(VGroup):
    """
    O que faz
    Desenha um arco marcador do ângulo de latitude em relação ao plano do horizonte, posicionando-o automaticamente no hemisfério correspondente.

    Como usar
        marcador = MarcadorLatitude(latitude=23.5, cor_arco=WHITE)
        cena.add(marcador)

    Parâmetros
        latitude (float, opcional): Latitude em graus que será marcada. Padrão é 0.
        cor_arco (Color, opcional): Cor do arco marcador. Padrão é WHITE.
        cor_texto (Color, opcional): Cor do texto do rótulo, se usado. Padrão é WHITE.
        espessura_arco (float, opcional): Largura do traço do arco. Padrão é 3.
        label_2d (bool, opcional): Exibe rótulo em 2D orientado para a câmera. Padrão é False.
        label_3d (bool, opcional): Exibe rótulo como objeto 3D na cena. Padrão é False.
        math_label_2d (str, opcional): Texto LaTeX personalizado para rótulo 2D. Padrão é None.
        math_label_3d (str, opcional): Texto LaTeX personalizado para rótulo 3D. Padrão é None.
        label_distance (float, opcional): Fator de distância do rótulo ao arco. Padrão é 1.
        tamanho_da_fonte (int, opcional): Tamanho da fonte do rótulo. Padrão é 20.
        arco (bool, opcional): Controla a adição do arco à cena. Padrão é True.
        **kwargs (dict, opcional): Argumentos adicionais herdados de VGroup.

    Observações
        Se latitude for zero, o marcador não desenha o arco, pois não há inclinação a ser marcada.
    """
    def __init__(self, latitude, cor_arco=WHITE, cor_texto=WHITE, espessura_arco=3, label_2d=False,label_3d=False, math_label_2d=None, math_label_3d=None, label_distance=1, tamanho_da_fonte=20, arco=True, **kwargs):

        super().__init__(**kwargs)
       
        Origem = ORIGIN
        
        # Obtém as posições dos dois pontos no espaço
        # Converte a latitude para radianos
        rad_latitude = np.radians(latitude)
        
        # Determina a posição do ponto no espaço 3D com base na latitude
        if latitude > 0:
            ponto1_pos = 2 * np.array([0, np.cos(rad_latitude), np.sin(rad_latitude)])
        else:
            ponto1_pos = 2 * np.array([0, -np.cos(rad_latitude), -np.sin(rad_latitude)])
        ponto2_pos = [ponto1_pos[0],ponto1_pos[1],0]

        # Se a latitude não for zero, cria o marcador


        
        if latitude!=0:

        
            # Obtém os vetores das duas linhas
            linha1_vector = ponto1_pos - Origem
            linha2_vector = ponto2_pos - Origem

            # Cálculo do ângulo entre os dois vetores usando o produto escalar
            cos_theta = np.dot(linha1_vector, linha2_vector) / (np.linalg.norm(linha1_vector) * np.linalg.norm(linha2_vector))
            angulo = np.arccos(cos_theta)  # Em radianos

            # Agora, criamos um arco entre os dois vetores no plano formado pelo produto vetorial
            # Vetor perpendicular ao plano formado pelos vetores linha1_vector e linha2_vector
            produto_vetorial = np.cross(linha1_vector, linha2_vector)
            produto_vetorial /= np.linalg.norm(produto_vetorial) 
            
            produto_vetorial2 = np.cross(linha1_vector, produto_vetorial)
            produto_vetorial2 /= np.linalg.norm(produto_vetorial2)

            num_pontos = 50
            pontos_arco = [
                    Origem + np.cos(angulo) * linha1_vector / np.linalg.norm(linha1_vector) * 0.5 + np.sin(angulo) * -produto_vetorial2 / np.linalg.norm(produto_vetorial2) * 0.5
                    for angulo in np.linspace(np.radians(0), angulo, num_pontos)
                ]
                
            arco2 = VMobject().set_stroke(width=espessura_arco,color=cor_arco)  # Chama a inicializaÃ§Ã£o de CirculoTangente
            arco2.set_points_as_corners(pontos_arco)
            
            if label_2d:
                texto = Tex(fr"{np.round(np.degrees(angulo),2)}^\circ", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
                texto.rotate(PI/2, X_AXIS).rotate(PI/2, Z_AXIS)
                self.add(texto)
            if label_3d:
                texto = Tex(fr"{np.round(np.degrees(angulo),2)}^\circ", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
                self.add(texto)
            if math_label_2d != None:
                texto = Tex(fr"{math_label_2d}", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
                texto.rotate(PI/2, X_AXIS).rotate(PI/2, Z_AXIS)
                self.add(texto)
            if math_label_3d != None:
                texto = Tex(fr"{math_label_3d}", font_size = tamanho_da_fonte).move_to(label_distance*1.7*arco.get_center_of_mass()).set_color(cor_texto)
                self.add(texto)
                
        

            # Adiciona as linhas e o arco ao VGroup
            
            self.add(arco2)

class MarcadorAnguloExterno(VGroup):
    """
    O que faz
    Constrói um arco externo entre dois pontos em coordenadas 3D sobre a esfera utilizando interpolação esférica e, opcionalmente, adiciona um rótulo.

    Como usar
        arco = MarcadorAnguloExterno(ponto1=obj1, ponto2=obj2, raio=2.2, cor=ORANGE, label="α")
        cena.add(arco)

    Parâmetros
        ponto1 (VMobject, opcional): Primeiro ponto sobre a esfera. Padrão é None.
        ponto2 (VMobject, opcional): Segundo ponto sobre a esfera. Padrão é None.
        raio (float, opcional): Raio efetivo para posicionar o arco. Padrão é 2.2.
        cor (Color, opcional): Cor do arco. Padrão é ORANGE.
        num_pontos (int, opcional): Número de amostras para o arco. Padrão é 50.
        label_distance (float, opcional): Fator de deslocamento do rótulo em relação ao centro de massa do arco. Padrão é 1.
        tamanho_fonte (int, opcional): Tamanho da fonte do rótulo. Padrão é 20.
        espessura (float, opcional): Largura do traço do arco. Padrão é 2.
        label (str, opcional): Texto do rótulo LaTeX a ser exibido. Padrão é None.

    Observações
        Se os pontos coincidirem, nada é desenhado. O rótulo é rotacionado para melhor legibilidade em cena 3D.
    """
    def __init__(self, ponto1:VMobject, ponto2:VMobject, raio=2.2, cor=ORANGE, num_pontos=50, label_distance=1, tamanho_fonte=20, espessura=2, label=None):
        super().__init__()
       
        inicio = ponto1.get_center()
        fim = ponto2.get_center()

        if np.allclose(inicio, fim):
            return
        
        angulo = np.arccos(np.dot(inicio, fim) / (np.linalg.norm(inicio) * np.linalg.norm(fim)))
        
        pontos_arco = []
        for t in range(num_pontos + 1):
            ponto_slerp = (
                np.sin((1 - t / num_pontos) * angulo) * inicio +
                np.sin((t / num_pontos) * angulo) * fim
            ) / np.sin(angulo) * raio / 2
            pontos_arco.append(ponto_slerp)
        
        arco = VMobject(color=cor, stroke_width=espessura)
        arco.set_points_as_corners(pontos_arco)
        
        if label is not None:
            texto = Tex(label, font_size=tamanho_fonte, color=cor).move_to(label_distance * 1.1 * arco.get_center_of_mass())
            texto.rotate(PI/2, X_AXIS).rotate(PI/2, Z_AXIS)
            self.add(texto)
    
        self.add(arco)

class SegmentoAstro(LinhaVMobject):
    """
    O que faz
    Desenha um segmento de linha do centro até um PontoAstro, útil para indicar a direção radial do astro na esfera.

    Como usar
        seg = SegmentoAstro(ponto=estrela, cor=WHITE, espessura=2)
        cena.add(seg)

    Parâmetros
        ponto (VMobject, opcional): Objeto que possui posição sobre a esfera, tipicamente PontoAstro. Padrão é None.
        cor (Color, opcional): Cor do segmento. Padrão é WHITE.
        espessura (float, opcional): Largura do traço. Padrão é 2.
        ORIGEM (VMobject ou np.ndarray, opcional): Objeto ou vetor indicando a origem do segmento. Padrão é ORIGIN.
        **kwargs (dict, opcional): Argumentos adicionais herdados de VMobject.

    Observações
        Se ORIGEM for um VMobject, usa seu centro como origem geométrica do segmento.
    """
    def __init__(self, ponto, cor=WHITE, espessura=2, ORIGEM = None, **kwargs):

        if ORIGEM == None:
            ORIGEM = ORIGIN
        else: 
            ORIGEM = ORIGEM.get_center()

        # Obtém a posição do ponto no espaço
        ponto_pos = ponto.get_center()
        super().__init__(ponto_final=ponto_pos, ponto_inicial=ORIGEM,espessura=espessura,cor=cor,**kwargs)



#Objetos 

class Moon(TexturedSurface):
    """
    O que faz
    Cria uma esfera texturizada representando a Lua, com controle de sombreamento e textura noturna.

    Como usar
        lua = Moon(radius=RENDER_MOON_RADIUS)
        cena.add(lua)

    Parâmetros
        radius (float, opcional): Raio da Lua. Padrão é RENDER_MOON_RADIUS.
        resolution (tuple, opcional): Resolução da esfera (u, v). Padrão é (101, 51).
        texture (str, opcional): URL ou caminho da textura diurna. Padrão é uma textura pública da Lua.
        dark_texture (str, opcional): URL ou caminho da textura noturna. Padrão é textura preta.
        shading (tuple, opcional): Parâmetros de sombreamento da superfície. Padrão é (0.25, 0.25, 1).
        **kwargs (dict, opcional): Argumentos adicionais herdados de TexturedSurface.

    Observações
        A rotação e posicionamento podem ser ajustados após a criação com os métodos padrão do Manim.
    """
    def __init__(
        self,
        radius=RENDER_MOON_RADIUS,
        resolution=(101, 51),
        texture="data_image/moon.jpg",
        dark_texture="data_image/black.jpg",
        shading=(0.25, 0.25, 1),
        **kwargs
    ):
        sphere = Sphere(radius=radius, resolution=resolution)
        super().__init__(sphere, texture, dark_texture, **kwargs)
        self.set_shading(*shading)

class Mars(TexturedSurface):
    """
    O que faz
    Cria uma esfera texturizada representando Marte, com controle de sombreamento e textura noturna.

    Como usar
        marte = Mars(radius=RENDER_MARS_RADIUS)
        cena.add(marte)

    Parâmetros
        radius (float, opcional): Raio de Marte. Padrão é RENDER_MARS_RADIUS.
        resolution (tuple, opcional): Resolução da esfera (u, v). Padrão é (101, 51).
        texture (str, opcional): URL ou caminho da textura diurna de Marte. Padrão é uma textura pública.
        dark_texture (str, opcional): URL ou caminho da textura noturna. Padrão é textura preta.
        shading (tuple, opcional): Parâmetros de sombreamento da superfície. Padrão é (0.25, 0.25, 1).
        **kwargs (dict, opcional): Argumentos adicionais herdados de TexturedSurface.

    Observações
        Sombreamento e opacidade podem ser ajustados para compor cenas realistas.
    """
    def __init__(
        self,
        radius=RENDER_MARS_RADIUS,
        resolution=(101, 51),
        texture="data_image/mars.jpg",
        dark_texture="data_image/black.jpg",
        shading=(0.25, 0.25, 1),
        **kwargs
    ):
        sphere = Sphere(radius=radius, resolution=resolution)
        super().__init__(sphere, texture, dark_texture, **kwargs)
        self.set_shading(*shading)

class Mercury(TexturedSurface):
    """
    O que faz
    Cria uma esfera texturizada representando Mercúrio, com controle de sombreamento e textura noturna.

    Como usar
        mercurio = Mercury(radius=RENDER_MERCURY_RADIUS)
        cena.add(mercurio)

    Parâmetros
        radius (float, opcional): Raio de Mercúrio. Padrão é RENDER_MERCURY_RADIUS.
        resolution (tuple, opcional): Resolução da esfera (u, v). Padrão é (101, 51).
        texture (str, opcional): Caminho da textura diurna. Padrão é "mercury.jpg".
        dark_texture (str, opcional): Caminho da textura noturna. Padrão é preto.
        shading (tuple, opcional): Parâmetros de sombreamento da superfície. Padrão é (0.25, 0.25, 1).
        **kwargs (dict, opcional): Argumentos adicionais herdados de TexturedSurface.
    """
    def __init__(
        self,
        radius=RENDER_MERCURY_RADIUS,
        resolution=(101, 51),
        texture="data_image/mercury.jpg",
        dark_texture="data_image/black.jpg",
        shading=(0.25, 0.25, 1),
        **kwargs
    ):
        sphere = Sphere(radius=radius, resolution=resolution)
        super().__init__(sphere, texture, dark_texture, **kwargs)
        self.set_shading(*shading)


class Venus(TexturedSurface):
    """
    O que faz
    Cria uma esfera texturizada representando Vênus, com controle de sombreamento e textura noturna.

    Como usar
        venus = Venus(radius=RENDER_VENUS_RADIUS)
        cena.add(venus)
    """
    def __init__(
        self,
        radius=RENDER_VENUS_RADIUS,
        resolution=(101, 51),
        texture="data_image/venus.jpg",
        dark_texture="data_image/black.jpg",
        shading=(0.25, 0.25, 1),
        **kwargs
    ):
        sphere = Sphere(radius=radius, resolution=resolution)
        super().__init__(sphere, texture, dark_texture, **kwargs)
        self.set_shading(*shading)


class Jupiter(TexturedSurface):
    """
    O que faz
    Cria uma esfera texturizada representando Júpiter, com controle de sombreamento e textura noturna.

    Como usar
        jupiter = Jupiter(radius=RENDER_JUPITER_RADIUS)
        cena.add(jupiter)
    """
    def __init__(
        self,
        radius=RENDER_JUPITER_RADIUS,
        resolution=(101, 51),
        texture="data_image/jupiter.jpg",
        dark_texture="data_image/black.jpg",
        shading=(0.25, 0.25, 1),
        **kwargs
    ):
        sphere = Sphere(radius=radius, resolution=resolution)
        super().__init__(sphere, texture, dark_texture, **kwargs)
        self.set_shading(*shading)


class Saturn(TexturedSurface):
    """
    O que faz
    Cria uma esfera texturizada representando Saturno, com controle de sombreamento e textura noturna.

    Como usar
        saturno = Saturn(radius=RENDER_SATURN_RADIUS)
        cena.add(saturno)
    """
    def __init__(
        self,
        radius=RENDER_SATURN_RADIUS,
        resolution=(101, 51),
        texture="data_image/saturn.jpg",
        dark_texture="data_image/black.jpg",
        shading=(0.25, 0.25, 1),
        **kwargs
    ):
        sphere = Sphere(radius=radius, resolution=resolution)
        super().__init__(sphere, texture, dark_texture, **kwargs)
        self.set_shading(*shading)


class Uranus(TexturedSurface):
    """
    O que faz
    Cria uma esfera texturizada representando Urano, com controle de sombreamento e textura noturna.

    Como usar
        urano = Uranus(radius=RENDER_URANUS_RADIUS)
        cena.add(urano)
    """
    def __init__(
        self,
        radius=RENDER_URANUS_RADIUS,
        resolution=(101, 51),
        texture="data_image/uranus.jpg",
        dark_texture="data_image/black.jpg",
        shading=(0.25, 0.25, 1),
        **kwargs
    ):
        sphere = Sphere(radius=radius, resolution=resolution)
        super().__init__(sphere, texture, dark_texture, **kwargs)
        self.set_shading(*shading)


class Neptune(TexturedSurface):
    """
    O que faz
    Cria uma esfera texturizada representando Netuno, com controle de sombreamento e textura noturna.

    Como usar
        netuno = Neptune(radius=RENDER_NEPTUNE_RADIUS)
        cena.add(netuno)
    """
    def __init__(
        self,
        radius=RENDER_NEPTUNE_RADIUS,
        resolution=(101, 51),
        texture="data_image/neptune.jpg",
        dark_texture="data_image/black.jpg",
        shading=(0.25, 0.25, 1),
        **kwargs
    ):
        sphere = Sphere(radius=radius, resolution=resolution)
        super().__init__(sphere, texture, dark_texture, **kwargs)
        self.set_shading(*shading)

class Sun(Group):
    """
    O que faz
    Agrupa uma superfície texturizada do Sol com halos de brilho próximos e amplos, oferecendo controles práticos de intensidade, tamanho e opacidade.

    Como usar
        sol = Sun(radius=RENDER_SUN_RADIUS, near_glow_ratio=0.2, big_glow_ratio=0.8)
        cena.add(sol)

    Parâmetros
        radius (float, opcional): Raio do Sol. Padrão é RENDER_SUN_RADIUS.
        texture (str, opcional): URL ou caminho da textura solar. Padrão é uma textura pública.
        near_glow_ratio (float, opcional): Razão do raio do halo interno em relação ao raio solar. Padrão é 0.0.
        near_glow_factor (float, opcional): Intensidade do brilho interno. Padrão é 5.
        big_glow_ratio (float, opcional): Razão do raio do halo externo em relação ao raio solar. Padrão é 0.7*RENDER_SUN_RADIUS.
        big_glow_factor (float, opcional): Intensidade do brilho externo. Padrão é 1.
        big_glow_opacity (float, opcional): Opacidade do brilho externo. Padrão é 0.35.
        shading (tuple, opcional): Parâmetros de sombreamento. Padrão é (0, 0, 0).
        edge (np.ndarray ou direção, opcional): Posição de ancoragem na tela via to_edge. Padrão é LEFT.
        **kwargs (dict, opcional): Argumentos adicionais herdados de TexturedSurface e Group.

    Observações
        Os métodos set_* permitem alterar dinamicamente brilho e texturas após a criação do objeto.
    """
    def __init__(
        self,
        radius=RENDER_SUN_RADIUS,
        texture="data_image/sun.jpg",
        near_glow_ratio=0.0,
        near_glow_factor=5,
        big_glow_ratio=0.7*RENDER_SUN_RADIUS,
        big_glow_factor=1,
        big_glow_opacity=0.35,
        shading=(0, 0, 0),
        edge=LEFT,
        **kwargs
    ):
        self.radius = radius
        self.texture = texture
        self.shading = shading
        self.edge = edge

        # Salva os valores dos parâmetros dos glows
        self.near_glow_ratio = near_glow_ratio
        self.near_glow_factor = near_glow_factor
        self.big_glow_ratio = big_glow_ratio
        self.big_glow_factor = big_glow_factor
        self.big_glow_opacity = big_glow_opacity

        # Criação dos componentes
        self.sun_surface = TexturedSurface(Sphere(radius=radius), texture, **kwargs)
        self.sun_surface.set_shading(*shading)
        self.sun_surface.to_edge(edge)

        self.near_glow = GlowDot(radius=near_glow_ratio * radius, glow_factor=near_glow_factor)
        self.near_glow.move_to(self.sun_surface)

        self.big_glow = GlowDot(radius=big_glow_ratio * radius, glow_factor=big_glow_factor, opacity=big_glow_opacity)
        self.big_glow.move_to(self.sun_surface)

        super().__init__(self.sun_surface, self.near_glow, self.big_glow)

    def set_near_glow_factor(self, factor):
        self.near_glow_factor = factor
        self.near_glow.set_glow_factor(factor)

    def set_near_glow_ratio(self, ratio):
        self.near_glow_ratio = ratio
        self.near_glow.set_radius(ratio * self.radius)

    def set_big_glow_factor(self, factor):
        self.big_glow_factor = factor
        self.big_glow.set_glow_factor(factor)

    def set_big_glow_ratio(self, ratio):
        self.big_glow_ratio = ratio
        self.big_glow.set_radius(ratio * self.radius)

    def set_big_glow_opacity(self, opacity):
        self.big_glow_opacity = opacity
        self.big_glow.set_opacity(opacity)

    def set_texture(self, texture):
        self.texture = texture
        self.sun_surface.set_texture(texture)

    def set_shading(self, shading):
        self.shading = shading
        self.sun_surface.set_shading(*shading)
        
class Earth(TexturedSurface):
    """
    O que faz
    Cria um globo terrestre texturizado com opção de nuvens animadas e texturas diurna e noturna.

    Como usar
        terra = Earth(radius=RENDER_EARTH_RADIUS, clouds=True)
        cena.add(terra)

    Parâmetros
        radius (float, opcional): Raio da Terra. Padrão é RENDER_EARTH_RADIUS.
        clouds (bool, opcional): Habilita camada de nuvens com leve rotação. Padrão é True.
        resolution (tuple, opcional): Resolução da esfera (u, v). Padrão é (101, 51).
        day_texture (str, opcional): URL ou caminho da textura diurna. Padrão é uma textura pública.
        night_texture (str, opcional): URL ou caminho da textura noturna. Padrão é uma textura pública.

    Observações
        O globo é rotacionado 90 graus em Z para alinhar mapas comuns com a cena.
    """
    def __init__(self, radius=RENDER_EARTH_RADIUS, clouds=True, resolution = (101,51),
                 day_texture="data_image/earth_day_high.jpg",
                 night_texture="data_image/earth_night.jpg"):
        sphere = Sphere(radius=radius,resolution=resolution)
        super().__init__(sphere, day_texture, night_texture)
        self.rotate(PI/2, Z_AXIS)
        if clouds:
            nuvem = Clouds(radius*1.02)
            nuvem.add_updater(lambda n,dt:n.rotate(0.1*dt,OUT))
            self.add(nuvem)

class Clouds(TexturedSurface):
    """
    O que faz
    Adiciona uma camada semitransparente de nuvens ao redor da Terra, com possibilidade de animação.

    Como usar
        nuvens = Clouds(radius=RENDER_EARTH_RADIUS*1.02)
        self.add(nuvens)

    Parâmetros
        radius (float, opcional): Raio da camada de nuvens. Padrão é RENDER_EARTH_RADIUS*1.02.
        day_texture (str, opcional): URL ou caminho da textura das nuvens. Padrão é um PNG transparente.
        night_texture (str, opcional): URL ou caminho da textura noturna das nuvens. Padrão é um PNG transparente.

    Observações
        A opacidade é reduzida para 0.2 e o sombreamento é neutralizado para efeito atmosférico suave.
    """
    def __init__(self, radius=RENDER_EARTH_RADIUS*1.02,
                 day_texture="data_image/clouds.png",
                 night_texture="data_image/clouds.png"):
        sphere = Sphere(radius=radius)
        super().__init__(sphere, day_texture, night_texture)
        self.set_shading(0,0,0).set_opacity(0.2)
        
        
        
#ADD

class PlanoRetangular(Surface):
    """
    O que faz
    Cria um plano retangular no plano XY como superfície paramétrica, útil como fundo, painel ou base para anotações.

    Como usar
        plano = PlanoRetangular(largura=10, altura=6, cor=BLUE, opacidade=0.5)
        cena.add(plano)

    Parâmetros
        largura (float, opcional): Extensão no eixo X do plano. Padrão é 10.
        altura (float, opcional): Extensão no eixo Y do plano. Padrão é 10.
        cor (Color, opcional): Cor de preenchimento do plano. Padrão é BLUE.
        opacidade (float, opcional): Transparência do preenchimento. Padrão é 1.
        resolucao (tuple, opcional): Resolução da malha paramétrica. Padrão é (20, 20).
        **kwargs (dict, opcional): Argumentos adicionais herdados de Surface.

    Observações
        A parametrização usa u em x e v em y, ambos no plano Z igual a zero.
    """
    def __init__(self, largura=10, altura=10, cor=BLUE, opacidade=1, resolucao=(20, 20), **kwargs):
        self.largura = largura
        self.altura = altura
        super().__init__(
            u_range=[-largura / 2, largura / 2],
            v_range=[-altura / 2, altura / 2],
            resolution=resolucao,
            color=cor,
            **kwargs
        )
        self.set_opacity(opacidade)

    def uv_func(self, u, v):
        return np.array([u, v, 0])
