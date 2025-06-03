from manimlib import *
import numpy as np

CELESTIAL_DEFAULT_RAIO = 6
SUN_DEFAULT_RADIUS = 10
MOON_DEFAULT_RADIUS = 1
DEFAULT_EARTH_RADIUS = 3
radius_orbit_earth = 50
radius_orbit_moon = 10

def angulo_diedro_de_vetores_cartesianos(A, B, V):
    """
        Encontra o ângulo diedro com vertice na posição V e os dois planos compostos pelo vértice V, origem e pontos A ou B
        (Angulo usado na trigonometria esférica como o ângulo interno de um triângulo esférico)

        Parâmetros:
            A (np.ndarray): Vetor de posição do primeiro ponto na esfera.
            B (np.ndarray): Vetor de posição do segundo ponto na esfera.
            V (np.ndarray): Vetor de posição do vértice onde o ângulo será calculado.

        Retorno:
            float: O ângulo esférico interno do triângulo esférico, em graus.
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

#ELEMENTOS BÁSICOS
class EsferaCeleste(Sphere):
    """
    Representa uma esfera celeste personalizada, estendendo a classe Sphere do Manim.

    Parâmetros:
        raio (float, opcional): Raio da esfera. Padrão é 2.
        cor (Color, opcional): Cor da esfera. Padrão é BLUE.
        cor_contorno (Color, opcional): Cor do contorno da esfera. Padrão é BLUE.
        resolucao (tuple, opcional): Resolução da esfera (quantidade de subdivisões). Padrão é (20, 20).
        opacidade (float, opcional): Opacidade do preenchimento da esfera. Padrão é 0.3.
        largura_contorno (float, opcional): Largura do contorno da esfera. Padrão é 3.
        opacidade_contorno (float, opcional): Opacidade do contorno da esfera. Padrão é 0.
        tabuleiro_checkered (tuple, opcional): Cores do padrão quadriculado da esfera. Padrão é (BLUE, BLUE).
        **kwargs: Argumentos adicionais para a classe Sphere.

    """
    def __init__(self, raio=CELESTIAL_DEFAULT_RAIO, cor=BLUE, resolucao=(80, 80), opacidade=0.3, **kwargs):
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

class SuperficieObservador(Surface):
    """
    Representa a superfície de observação como um disco plano no plano XY, 
    utilizando a classe Surface do Manim.

    Parâmetros:
        raio (float, opcional): Raio do disco que representa a superfície de observação. Padrão é 2.
        cor_preenchimento (Color, opcional): Cor do preenchimento da superfície. Padrão é GREEN.
        opacidade (float, opcional): Opacidade do preenchimento. Padrão é 1 (totalmente visível).
        resolucao (tuple, opcional): Resolução da malha da superfície. Padrão é (20, 20).
        **kwargs: Argumentos adicionais para a classe Surface.

    """
    def __init__(self, raio=CELESTIAL_DEFAULT_RAIO, cor=GREEN,resolucao=(50,50),opacidade=1, **kwargs):
        self.raio = raio
        super().__init__(
            u_range=[0, TAU],
            v_range=[0, 1],
            color=cor,
            resolution=resolucao,
            **kwargs
        )
        self.set_opacity(opacidade)

    def uv_func(self, u, v):
        return np.array([
            self.raio * np.cos(u) * v,
            self.raio * np.sin(u) * v,
            0
        ])   
        
class PontoAstro(Sphere):
    """
    Representa um ponto na esfera celeste a partir de coordenadas esféricas (altura e azimute),
    herdando da classe Dot3D do Manim.

    OBS.: A direção Norte é o vetor espacial [0,1,0]

    Parâmetros:
        altura (float): Altura do ponto (equivalente à latitude esférica, em graus, de -90 a +90).
        azimute (float): Azimute do ponto (equivalente à longitude esférica, em graus, de 0 a 360 partindo do Norte pra direção Leste).
        tamanho (float, opcional): Raio do ponto 3D. Padrão é DEFAULT_DOT_RADIUS.
        cor (Color, opcional): Cor do ponto. Padrão é WHITE.
        raio (float, opcional): Raio da esfera na qual o ponto será posicionado. Padrão é 2.
        **kwargs: Argumentos adicionais para a classe Dot3D.

    """
    def __init__(self, altura, azimute, tamanho=DEFAULT_DOT_RADIUS*0.6, cor=WHITE, raio=CELESTIAL_DEFAULT_RAIO, center = None, **kwargs):
        
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
        if center == None:
            self.move_to(np.array([x, y, z]))
        else:
            self.move_to(np.array([x, y, z])+center)
            

class P(PontoAstro):
    def __init__(self, altura, azimute, tamanho=DEFAULT_DOT_RADIUS*0.6, cor=WHITE, raio=CELESTIAL_DEFAULT_RAIO, center =None,**kwargs):
        super().__init__(altura, azimute, tamanho=tamanho, cor=cor, raio=raio,center=center, **kwargs)

class PontoAstroEquatorial(PontoAstro):
    """
    Representa um ponto na esfera celeste a partir de coordenadas esféricas (altura e azimute),
    herdando da classe Dot3D do Manim.

    OBS.: A direção Norte é o vetor espacial [0,1,0]

    Parâmetros:
        altura (float): Altura do ponto (equivalente à latitude esférica, em graus, de -90 a +90).
        azimute (float): Azimute do ponto (equivalente à longitude esférica, em graus, de 0 a 360 partindo do Norte pra direção Leste).
        tamanho (float, opcional): Raio do ponto 3D. Padrão é DEFAULT_DOT_RADIUS.
        cor (Color, opcional): Cor do ponto. Padrão é WHITE.
        raio (float, opcional): Raio da esfera na qual o ponto será posicionado. Padrão é 2.
        **kwargs: Argumentos adicionais para a classe Dot3D.

    """
    def __init__(self, declinacao, ascencao_reta_graus, latitude, TSL_graus=0,  tamanho=DEFAULT_DOT_RADIUS*0.6, cor=WHITE, raio=CELESTIAL_DEFAULT_RAIO, **kwargs):
        
        # Armazena os atributos do ponto
        self.declinacao = declinacao  # Equivalente à latitude esférica
        self.ascencao_reta = ascencao_reta_graus  # Equivalente à longitude esférica
        self.cor = cor
        self.raio = raio
        self.tamanho = tamanho


        # Converte a altura e a longitude para radianos
         # Ângulo polar (90° - altura para alinhar com a convenção esférica)
        phi = self.ascencao_reta + TSL_graus + 180  # Ângulo azimutal
        
        # Converte as coordenadas esféricas para coordenadas cartesianas
        
        # Chama o construtor da classe pai (Dot3D) para criar o ponto na posição calculada
        super().__init__(self.declinacao,phi,tamanho=self.tamanho, cor=self.cor, raio=raio,**kwargs)
        
        self.rotate_about_origin((latitude - 90) * DEGREES, X_AXIS)

class PontoVernal(PontoAstroEquatorial):
    def __init__(self, latitude, TSL_graus=0,  tamanho=DEFAULT_DOT_RADIUS*0.6, cor=PINK, raio=CELESTIAL_DEFAULT_RAIO, **kwargs):
        super().__init__(0,0,latitude=latitude,TSL_graus=TSL_graus,tamanho=tamanho, cor=cor, raio=raio,**kwargs)

class GrandeArco(VMobject):
    """
    Representa um grande arco esférico entre dois pontos em uma esfera,
    utilizando interpolação esférica (Slerp) para gerar uma curva suave.

    Parâmetros:
        ponto1 (Dot3D): Primeiro ponto do arco.
        ponto2 (Dot3D): Segundo ponto do arco.
        raio (float, opcional): Raio da esfera sobre a qual o arco será traçado. Padrão é 2.
        cor (Color, opcional): Cor do arco. Padrão é WHITE.
        num_pontos (int, opcional): Número de pontos para interpolação do arco. Padrão é 50.
        espessura (float, opcional): Espessura da linha do arco. Padrão é 2.
    
    """
    def __init__(self, ponto1, ponto2, cor=WHITE, num_pontos=50, center = None, espessura=2):
        # Inicializa o objeto como um VGroup (grupo de vetores gráficos do Manim)
        super().__init__()
        if center==None:
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
        
class GrandeCirculo(ParametricCurve):
    def __init__(self, vetor_normal, raio=CELESTIAL_DEFAULT_RAIO, cor=YELLOW, espessura=4):
        """
        Desenha um círculo máximo (grande círculo) em uma esfera dado um vetor normal ao plano do círculo.
        
        Parâmetros:
            vetor_normal (np.array): Vetor unitário 3D normal ao plano do grande círculo.
            raio (float): Raio da esfera (padrão: 2).
            cor (Color): Cor do grande círculo (padrão: BLUE).
            espessura (float): Espessura do círculo (padrão: 5).
        """
        vetor_normal = vetor_normal / np.linalg.norm(vetor_normal)
        raio_corrigido = raio + 0.02
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

        # Adiciona a função à instância da classe

class Paralelo(VMobject):
    def __init__(self, altitude, latitude=0, raio=CELESTIAL_DEFAULT_RAIO, num_pontos=400, cor=BLUE_C):
        """
        Desenha um círculo paralelo (de latitude) em uma esfera com o raio fornecido.

        Parâmetros:
            raio (float): Raio da esfera.
            altitude (float): altitude fixa (em graus) onde o círculo será desenhado.
            latitude (float): Latitude do lugar (Padrão = 0)
            num_pontos (int): Número de pontos para desenhar o círculo.
            cor (color): Cor do círculo.

        Retorna:
            VMobject: Um círculo na latitude fornecida.
        """
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
        self.set_points_as_corners(pontos_círculo).set_color(cor)
        
        # Adicionar o círculo ao VGroup
        if latitude != 0:
            self.rotate_about_origin((latitude - 90) * DEGREES, X_AXIS)

class ArcoParalelo(Paralelo):
    def __init__(self, altitude, angulo_inicial, angulo_final, raio=CELESTIAL_DEFAULT_RAIO, num_pontos=50, cor=WHITE):
        """
        Desenha um arco de círculo paralelo (de latitude) em uma esfera com o raio fornecido.

        Parâmetros:
            altitude (float): Latitude fixa (em graus) onde o arco será desenhado.
            angulo_inicio (float): Azimute inicial do arco (em graus).
            angulo_fim (float): Azimute final do arco (em graus).
            raio (float): Raio da esfera.
            num_pontos (int): Número de pontos para desenhar o arco.
            cor (color): Cor do arco.
        """
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
        self.set_points_as_corners(pontos_arco).set_color(cor)
        
class MeridianoLocal(Arc):
    def __init__(self, raio=CELESTIAL_DEFAULT_RAIO, espessura=3, cor=ORANGE, opacidade=1, **kwargs):
        """
        Cria um meridiano local como um arco de círculo maior.

        raio: Raio do arco (padrão: 2)
        espessura: Largura da linha do arco (padrão: 3)
        cor: Cor do arco (padrão: WHITE)
        opacidade: Transparência do arco (padrão: 1)
        """
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
 
 
#DEPENDENTE DO OBSERVADOR
class EixoPolar(Line3D):
    """
    Cria um eixo polar em 3D a partir de uma latitude específica, representado por pequenos segmentos de linha.

    Parâmetros:
        latitude_graus (float): Latitude em graus.
        comprimento (float): Comprimento total do eixo (padrão: 4).
        phi (float): Ângulo de orientação azimutal do eixo (padrão: PI/2).
        espessura (float): Espessura dos segmentos de linha (padrão: 0.02).
        cor (Color): Cor dos segmentos de linha (padrão: RED).
        granularidade_3d (int): Número de segmentos para aumentar a suavidade (padrão: 20).
    """
    def __init__(self, latitude_graus, comprimento=2*CELESTIAL_DEFAULT_RAIO, phi=PI / 2, espessura=0.03, cor=RED):

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
    def __init__(self, latitude, cor=YELLOW_D, raio=CELESTIAL_DEFAULT_RAIO):
        """
        Cria um equador em uma esfera, dado um valor de latitude.
        
        Parâmetros:
            latitude (float): Latitude onde o equador será desenhado em graus.
            cor (Color): Cor do círculo (padrão: YELLOW_D).
            raio (float): Raio da esfera (padrão: 2).
        """
        # Converte a latitude para radianos
        latitude_rad = np.radians(latitude)
        
        
        
        # Define o vetor normal ao plano do equador
        vetor_normal = np.array([0, np.cos(latitude_rad), np.sin(latitude_rad)])
        
        # Chama o construtor da classe pai para criar o grande círculo
        super().__init__(vetor_normal=vetor_normal, raio=raio, cor=cor)

class Grade(VGroup):
    def __init__(self, latitude=90, raio=CELESTIAL_DEFAULT_RAIO, cor_dec=GREY_A,cor_ar=GREY_A, opacidade=0.8,espessura=2, n_ar=24,**kwargs):
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

class GradeMesh(SurfaceMesh):
    def __init__(self, latitude=90,raio=CELESTIAL_DEFAULT_RAIO, cor=BLUE_C, n_dec=19,n_ra=25, opacidade=0.5, espessura=2,**kwargs):
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



#ELEMENTOS EXTRA

class CirculoTangente(Polygon):
    def __init__(self, ponto: VMobject, raio_circulo=0.2, cor=WHITE, espessura=3):
        """
        Cria um círculo tangente à superfície de uma esfera no ponto 'p', com o ponto 'p' como centro.
        
        Parâmetros:
            ponto (Dot3D): O ponto na superfície da esfera onde o círculo é desenhado.
            raio_circulo (float): Raio do círculo tangente (padrão: 0.2).
            cor (Color): Cor do círculo (padrão: WHITE).
            espessura (float): Espessura do círculo (padrão: 1).
        """
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

class AnguloEsferico(VGroup):
    def __init__(self, p: VMobject, p1: VMobject, p2: VMobject, raio_circulo=0.3, cor=WHITE, math_label = None, cor_do_texto=WHITE, label_distance=1,tamanho_da_fonte=30, num_pontos=30, center=None,espessura=2):
        """
        Desenha um arco de cÃ­rculo tangente Ã  esfera no ponto 'p'.
        
        ParÃ¢metros:
            p (Dot3D): O ponto na superfÃ­cie da esfera.
            p1 (Dot3D): O primeiro ponto usado para calcular o Ã¢ngulo esfÃ©rico.
            p2 (Dot3D): O segundo ponto usado para calcular o Ã¢ngulo esfÃ©rico.
            raio_circulo (float): Raio do cÃ­rculo tangente (padrÃ£o: 0.3).
            cor (Color): Cor do arco (padrÃ£o: WHITE).
            num_pontos (int): NÃºmero de pontos para desenhar o arco (padrÃ£o: 30).
        
        """
        
        if center == None:
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
        normal_vector = p_pos / np.linalg.norm(p_pos)  # Vetor normal ao plano tangente

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
            center + 1.01*p_pos + raio_circulo * (np.cos(angulo) * -vetor_tangente2 + np.sin(angulo) * -vetor_tangente1) - normal_vector / 30
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
    
    def valor_angulo(self) -> float:  # Agora Ã© acessado como um atributo

        return angulo_diedro_de_vetores_cartesianos(self.p1.get_center(), self.p2.get_center(), self.p_pos)
   
    def atualizar_angulo(self, p: VMobject, p1: VMobject, p2: VMobject):
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
    def __init__(self, lat_ini, lat_fim, raio=CELESTIAL_DEFAULT_RAIO+0.02, cor=YELLOW, opacidade=0.5, resolucao=(40, 40), **kwargs):
        """
        Cria uma faixa de uma esfera delimitada pelas latitudes lat_ini e lat_fim.

        lat_ini: Latitude inicial em graus
        lat_fim: Latitude final em graus
        raio: Raio da esfera (padrão: 1)
        cor: Cor da faixa (padrão: YELLOW)
        opacidade: Transparência da faixa (padrão: 1)
        resolucao: Número de divisões na malha (padrão: 30x60)
        """
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
    def __init__(self, lon_ini, lon_fim, raio=CELESTIAL_DEFAULT_RAIO+0.02, cor=GREEN, opacidade=0.5, resolucao=(40, 40), **kwargs):
        
        """
        Cria uma faixa de uma esfera delimitada pelas longitudes lon_ini e lon_fim.

        lon_ini: Longitude inicial em graus
        lon_fim: Longitude final em graus
        raio: Raio da esfera (padrão: 2)
        cor: Cor da faixa (padrão: GREEN)
        opacidade: Transparência da faixa (padrão: 0.2)
        resolucao: Número de divisões na malha (padrão: 20x10)
        """
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
            def __init__(self, ponto_inicial, ponto_final, cor=WHITE, espessura=2, **kwargs):
                super().__init__(**kwargs)
                self.set_points_as_corners([ponto_inicial,(ponto_inicial+ponto_final)/2, ponto_final])
                self.set_color(cor)
                self.set_stroke(width=espessura)

#MARCADORES

class MarcadorAngulo(VGroup):
    def __init__(self, ponto1:VMobject, ponto2:VMobject, Origem = None , barra=True, cor_linha=WHITE, espessura_linha=2, cor_arco=WHITE, cor_texto=WHITE, espessura_arco=3, label_2d=False,label_3d=False, math_label_2d=None, math_label_3d=None, label_distance=1, tamanho_da_fonte=20, arco=True, **kwargs):
        """
        Classe que cria um arco marcador de ângulo entre duas linhas radiais que conectam o centro aos pontos.
        
        Parâmetros:
            ponto1 (PontoAstro): O primeiro ponto no espaço.
            ponto2 (PontoAstro): O segundo ponto no espaço.
            barra (bool): Se True, adiciona as linhas radiais.
            cor_linha (Color): Cor das linhas radiais.
            espessura_linha (float): Espessura das linhas radiais.
            cor_arco (Color): Cor do arco do ângulo.
            espessura_arco (float): Espessura do arco.
            label_2d (bool): Se True, exibe um rótulo 2D com o valor do ângulo.
            label_3d (bool): Se True, exibe um rótulo 3D com o valor do ângulo.
            math_label_2d (str): Texto matemático opcional para exibir como rótulo 2D.
            math_label_3d (str): Texto matemático opcional para exibir como rótulo 3D.
            label_distance (float): Distância do rótulo ao arco.
            tamanho_da_fonte (int): Tamanho da fonte do rótulo.
            arco (bool): Se True, adiciona o arco do ângulo.
            Origem (np.array): Ponto de origem para as linhas radiais.
        """
        super().__init__(**kwargs)
        if Origem == None:
            Origem = ORIGIN
        else: 
            Origem = Origem.get_center()
        # Obtém as posições dos dois pontos no espaço
        ponto1_pos = ponto1.get_center()
        ponto2_pos = ponto2.get_center()
        
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
            self.add(segmentos_linha1, segmentos_linha2)
        if arco:
            self.add(arco2)

class MarcadorAltura(VGroup):
    def __init__(self, ponto:VMobject, Origem = None , barra=True, cor_linha=WHITE, espessura_linha=2, cor_arco=WHITE, cor_texto=WHITE, espessura_arco=3, label_2d=False,label_3d=False, math_label_2d=None, math_label_3d=None, label_distance=1, tamanho_da_fonte=20, arco=True, **kwargs):
        """
        Classe que cria um arco marcador de ângulo entre duas linhas radiais que conectam o centro aos pontos.
        
        Parâmetros:
            ponto1 (PontoAstro): O primeiro ponto no espaço.
            ponto2 (PontoAstro): O segundo ponto no espaço.
            barra (bool): Se True, adiciona as linhas radiais.
            cor_linha (Color): Cor das linhas radiais.
            espessura_linha (float): Espessura das linhas radiais.
            cor_arco (Color): Cor do arco do ângulo.
            espessura_arco (float): Espessura do arco.
            label_2d (bool): Se True, exibe um rótulo 2D com o valor do ângulo.
            label_3d (bool): Se True, exibe um rótulo 3D com o valor do ângulo.
            math_label_2d (str): Texto matemático opcional para exibir como rótulo 2D.
            math_label_3d (str): Texto matemático opcional para exibir como rótulo 3D.
            label_distance (float): Distância do rótulo ao arco.
            tamanho_da_fonte (int): Tamanho da fonte do rótulo.
            arco (bool): Se True, adiciona o arco do ângulo.
            Origem (np.array): Ponto de origem para as linhas radiais.
        """
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
    Classe que cria um arco entre dois pontos no plano 2D e, opcionalmente, uma etiqueta com o valor do ângulo.

    Parâmetros:
        ponto1 (tuple): Coordenadas (x, y) do primeiro ponto no plano. Usado para definir o início do arco. 
        ponto2 (tuple): Coordenadas (x, y) do segundo ponto no plano. Usado para definir o final do arco. 
        raio (float): Raio do arco que conecta os dois pontos. Controla o tamanho do arco. 
        cor (cor): Cor do arco. Também define a cor da etiqueta se fornecida. 
        num_pontos (int): Número de pontos usados para aproximar a curva do arco. 
        label_distance (float): Distância entre o centro do arco e a etiqueta de texto.  
        tamanho_fonte (int): Tamanho da fonte da etiqueta de texto. 
        espessura (float): Espessura da linha do arco. Controla a espessura do arco gerado. 
        label (str, opcional): Texto a ser exibido como rótulo no arco. Se fornecido, uma etiqueta será adicionada ao arco. 
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
    def __init__(self, ponto, cor=WHITE, espessura=2, ORIGEM = None, **kwargs):
        """
        Classe que cria uma linha até o ponto 
        Parâmetros:
            ponto (PontoAstro): O ponto no espaço.
            cor_linha (cor): A cor das linhas radiais.
            espessura_linha (float): A espessura das linhas radiais.
        """
        if ORIGEM == None:
            ORIGEM = ORIGIN
        else: 
            ORIGEM = ORIGEM.get_center()

        # Obtém a posição do ponto no espaço
        ponto_pos = ponto.get_center()
        super().__init__(ponto_final=ponto_pos, ponto_inicial=ORIGEM,espessura=espessura,cor=cor,**kwargs)



#Objetos 

class Moon(TexturedSurface):
    def __init__(
        self,
        radius=MOON_DEFAULT_RADIUS,
        resolution=(101, 51),
        texture="https://upload.wikimedia.org/wikipedia/commons/thumb/7/74/Moon_texture.jpg/2560px-Moon_texture.jpg",
        dark_texture="https://upload.wikimedia.org/wikipedia/commons/thumb/5/50/Black_Wallpaper.jpg/1200px-Black_Wallpaper.jpg",
        shading=(0.25, 0.25, 1),
        **kwargs
    ):
        sphere = Sphere(radius=radius, resolution=resolution)
        super().__init__(sphere, texture, dark_texture, **kwargs)
        self.set_shading(*shading)

class Sun(Group):
    def __init__(
        self,
        radius=SUN_DEFAULT_RADIUS,
        texture="https://upload.wikimedia.org/wikipedia/commons/thumb/9/99/Map_of_the_full_sun.jpg/1280px-Map_of_the_full_sun.jpg",
        near_glow_ratio=0.0,
        near_glow_factor=5,
        big_glow_ratio=0.7*SUN_DEFAULT_RADIUS,
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
    def __init__(self, radius=DEFAULT_EARTH_RADIUS,
                 day_texture="https://upload.wikimedia.org/wikipedia/commons/thumb/4/4d/Whole_world_-_land_and_oceans.jpg/1280px-Whole_world_-_land_and_oceans.jpg",
                 night_texture="https://upload.wikimedia.org/wikipedia/commons/thumb/b/ba/The_earth_at_night.jpg/1280px-The_earth_at_night.jpg"):
        sphere = Sphere(radius=radius)
        super().__init__(sphere, day_texture, night_texture)
class Clouds(TexturedSurface):
    def __init__(self, radius=DEFAULT_EARTH_RADIUS+0.08,
                 day_texture="https://www.nicepng.com/png/full/120-1200066_earth-clouds-png-banner-library-library-earth-clouds.png",
                 night_texture="https://www.nicepng.com/png/full/120-1200066_earth-clouds-png-banner-library-library-earth-clouds.png"):
        sphere = Sphere(radius=radius).set_shading(0,0,0)
        super().__init__(sphere, day_texture, night_texture)
        
        
        
#ADD

class PlanoRetangular(Surface):
    """
    Representa um plano retangular no plano XY, utilizando a classe Surface do ManimGL.

    Parâmetros:
        largura (float, opcional): Tamanho no eixo X. Padrão é 4.
        altura (float, opcional): Tamanho no eixo Y. Padrão é 4.
        cor (Color, opcional): Cor do plano. Padrão é BLUE.
        opacidade (float, opcional): Opacidade do plano. Padrão é 1.
        resolucao (tuple, opcional): Resolução da malha da superfície. Padrão é (20, 20).
        **kwargs: Argumentos adicionais para a classe Surface.
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
