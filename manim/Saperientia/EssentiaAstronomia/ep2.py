from manimlib import *
from manimlib.Saperientia.astro_structures import * 
from manimlib.Saperientia.stellarium import * 
from manimlib.Saperientia.Variables import * 
from manimlib.Saperientia.camera import * 
import numpy as np
from scipy.spatial.transform import Rotation as R

class EclipseLunar(Scene):
    def construct(self):
        self.frame.reorient(80,80,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1.5)
        stars_data = extract_star_data()
        stars = Group()
        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            stars.add(star)
        self.add(stars)

        self.frame.reorient(90,45,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*2)
        terra = Earth().rotate(180*DEGREES,axis=OUT)
        sun = Sun(big_glow_ratio=0.2).scale(0.1).move_to([RENDER_EARTH_RADIUS*100,0,0])
        self.camera.light_source.move_to(sun.get_center())
        moon = Moon().move_to([-RENDER_EARTH_RADIUS*4,0,0]).rotate(180*DEGREES,axis=OUT)
        self.add(terra,sun,moon)
        self.wait(5)
        self.play(self.frame.animate.reorient(270,100,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3),run_time=7)
        
class CruzeiroSulMovimento(Scene):
    def construct(self):
        stars_data = extract_star_data()
        cruzeiro = Group()
        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            cruzeiro.add(star)
        self.add(cruzeiro)

        earth = Earth(radius=RENDER_EARTH_RADIUS, clouds=False, resolution=(400, 400))
        self.add(earth)

        # Observador: ponto na superfície
        lat = ValueTracker(25)  # Começa em São Paulo
        lon = 180
        R = RENDER_EARTH_RADIUS * 1.0001

        def camera_updater(mob):
            point = _latlon_to_xyz(lat.get_value(), lon, R)
            reorientar_camera(mob, point=point, theta_deg=90+lon, phi_deg=-lat.get_value()-20+180, gamma_deg=0, fovy=60)
        self.camera.frame.add_updater(camera_updater)

        # Anima o observador indo para o sul
        self.play(lat.animate.increment_value(-40), run_time=10, rate_func=linear)

        self.camera.frame.clear_updaters()
        self.play(AnimarCamera(self.camera.frame, np.array(self.camera.frame.get_implied_camera_location()) * 3 + np.array([0.01,0,0.1]), phi_deg=155, theta_deg=30, fovy=85))

        obs = PontoAstro(0, 235, 0.2 * ASTRO_SIZE, PURPLE,raio=RENDER_EARTH_RADIUS*1.01)
        ponto = TrueDot(obs.get_center(), radius = 0.2 * ASTRO_SIZE)
        self.add(ponto)
        lat_obs = ValueTracker(0)
        ponto.add_updater(lambda m: m.become(TrueDot(PontoAstro(lat_obs.get_value(), 235, 0.2 * ASTRO_SIZE, PURPLE,raio=RENDER_EARTH_RADIUS*1.01).get_center(), radius = 0.2 * ASTRO_SIZE, color = PURPLE, opacity = 0.7)))
        self.play(lat_obs.animate.increment_value(-30), run_time=10, rate_func=linear)
        
class TerraGirando(Scene):
    def construct(self):
        self.frame.reorient(90,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1.5)
        self.frame.set_field_of_view(-80)
        stars_data = extract_star_data()
        cruzeiro = Group()
        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            cruzeiro.add(star)
        self.add(cruzeiro)

        earth = Earth(radius=RENDER_EARTH_RADIUS, clouds=True)
        self.add(earth)
        earth.add_updater(lambda m,dt:m.rotate(dt*0.3,axis=OUT,about_point=ORIGIN))
        self.wait(5)
        self.play(self.frame.animate.reorient(270,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1.5),run_time=10,rate_func=linear)
        
class TerrasEntreLua(Scene):
    def construct(self):
        stars_data = extract_star_data()
        cruzeiro = Group()
        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            cruzeiro.add(star)
        self.add(cruzeiro)

        earth = Earth(radius=RENDER_EARTH_RADIUS, clouds=False)
        self.add(earth)
        lua = Moon()
        self.add(lua)
        lua.move_to(RENDER_EARTH_RADIUS*RIGHT*15.5)
        self.frame.reorient(0,90,0,lua.get_center()/2,RENDER_EARTH_RADIUS*12)
        self.camera.light_source.move_to(self.frame.get_implied_camera_location())
        self.wait(5)
        for i in range(1,8):
            earth.copy()
            self.play(ShowCreation(earth.copy().move_to(RENDER_EARTH_RADIUS*RIGHT*2*i).set_opacity(0.4)),run_time=0.3)
        self.wait(5)

class TerrasEntreLua30(Scene):
    def construct(self):
        stars_data = extract_star_data()
        cruzeiro = Group()
        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            cruzeiro.add(star)
        self.add(cruzeiro)

        earth = Earth(radius=RENDER_EARTH_RADIUS, clouds=False)
        self.add(earth)
        lua = Moon()
        self.add(lua)
        lua.move_to(RENDER_EARTH_RADIUS*RIGHT*61.5)
        self.frame.reorient(0,90,0,lua.get_center()/2,RENDER_EARTH_RADIUS*40)
        self.camera.light_source.move_to(self.frame.get_implied_camera_location())
        self.wait(5)
        for i in range(1,31):
            earth.copy()
            self.play(ShowCreation(earth.copy().move_to(RENDER_EARTH_RADIUS*RIGHT*2*i).set_opacity(0.4)),run_time=0.1)
        self.wait(5)

class FlatEarth(Surface):
    def __init__(
        self,
        radius: float = 1,
        u_range = (0, TAU),
        v_range = (1, 0),
        resolution = (2, 100),
        **kwargs
    ):
        super().__init__(
            u_range=u_range,
            v_range=v_range,
            resolution=resolution,
            **kwargs,
        )
        self.scale(radius)

    def uv_func(self, u: float, v: float) -> np.ndarray:
        return np.array([
            v * math.cos(u),
            v * math.sin(u),
            0
        ])

class TerraPlana(Scene):
    def construct(self):
        self.frame.reorient(70,45,0,ORIGIN,7)

        stars_data = extract_star_data()
        cruzeiro = Group()
        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            cruzeiro.add(star)
        self.add(cruzeiro.rotate_about_origin(-30*DEGREES))
        
        
        disk = FlatEarth(radius=3, resolution=(102, 100))
        sphere = Sphere(radius=3, resolution=(102,100))
        
        day_texture = "https://upload.wikimedia.org/wikipedia/commons/thumb/4/4d/Whole_world_-_land_and_oceans.jpg/1280px-Whole_world_-_land_and_oceans.jpg"
        night_texture = "https://upload.wikimedia.org/wikipedia/commons/thumb/b/ba/The_earth_at_night.jpg/1280px-The_earth_at_night.jpg"

        surfaces = [
            TexturedSurface(disk, day_texture, day_texture),
            TexturedSurface(sphere, day_texture, night_texture)
        ]

        ponto1 = PontoAstro(-30, -30, raio=3, tamanho=0.1,cor=PURPLE).set_z_index(5)
        ponto2 = PontoAstro(-30, -110, raio=3, tamanho=0.1,cor=PURPLE).set_z_index(5)

        seta1_final = Arrow(ponto1.get_center(), ponto1.get_center()*0.8 + 2.2*IN, thickness=8,buff=0,tip_angle=PI/3,color=YELLOW)
        seta2_final = Arrow(ponto2.get_center(), ponto2.get_center()*0.8 + 2.2*IN, thickness=8,buff=0,tip_angle=PI/3,color=YELLOW)
        
        c1 = ponto1.get_center()
        c2 = ponto2.get_center()


        ponto1.move_to(np.array([c1[0],c1[1], 0])/1.3)
        ponto2.move_to(np.array([c2[0],c2[1], 0])/1.3)

        seta1 = Arrow(ponto1.get_center(), ponto1.get_center()*1.4 + OUT*0.2, thickness=10,buff=0,tip_angle=PI/3,color=YELLOW)
        seta2 = Arrow(ponto2.get_center(), ponto2.get_center()*1.4 + OUT*0.2, thickness=10,buff=0,tip_angle=PI/3,color=YELLOW)



        self.play(FadeIn(surfaces[0]),)
        self.wait(3)
        self.add(ponto1, ponto2)
        self.play(FadeIn(ponto1),FadeIn(ponto2))
        self.wait(1)
        self.play(FadeIn(seta1),FadeIn(seta2))

        self.wait(3)
        self.play(self.frame.animate.reorient(250,60,0,ORIGIN,7),run_time=2)
        self.wait(1)

        self.play(
            Transform(surfaces[0], surfaces[1]),
            ponto1.animate.move_to(PontoAstro(-30, -30, raio=3, tamanho=0.1).get_center()),
            ponto2.animate.move_to(PontoAstro(-30, -110, raio=3, tamanho=0.1).get_center()),
            Transform(seta1, seta1_final),
            Transform(seta2, seta2_final),
            self.frame.animate.reorient(240,140,0,ORIGIN,8),
            run_time=3
        )
        self.wait(1)
        self.play(self.frame.animate.reorient(250-180,20,0,ORIGIN+(c1+c2),20),run_time=2)
        self.wait(4)

class FasesLua(Scene):
    def construct(self):
        stars_data = extract_star_data()
        stars = Group()
        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            stars.add(star)
        self.add(stars)

        self.frame.reorient(300,80,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3)
        terra = Earth().rotate(180*DEGREES,axis=OUT)
        sun = Sun(big_glow_ratio=0.1).scale(0.1).move_to([RENDER_EARTH_RADIUS*100,0,0])
        self.camera.light_source.move_to(sun.get_center())
        moon = Moon().move_to([-RENDER_EARTH_RADIUS*6,0,0]).rotate(180*DEGREES,axis=OUT)
        self.add(terra,sun,moon)
        self.wait(5)
        self.play(self.frame.animate.reorient(420,80,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3),run_time=2)
        self.wait(0.2)
        self.play(self.frame.animate.reorient(420,80,0,moon.get_center(),CELESTIAL_SPHERE_RADIUS*1),run_time=2)
        self.wait(0.2)
        self.play(self.frame.animate.reorient(420,80,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3),run_time=2)
        self.wait(1)
        angle = ValueTracker(0)
        
        moon.add_updater(lambda m:m.move_to([-np.cos(angle.get_value())*RENDER_EARTH_RADIUS*6,-np.sin(angle.get_value())*RENDER_EARTH_RADIUS*6,0]))
        moon.add_updater(lambda m:self.add(moon))
        self.play(angle.animate.increment_value(PI),Rotate(moon,PI,axis=OUT),run_time=2)        
        self.play(self.frame.animate.reorient(300,70,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3),run_time=3)
        self.wait(3)
        self.play(self.frame.animate.reorient(360,60,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*2.5),run_time=3)
        self.play(angle.animate.increment_value(PI),Rotate(moon,PI,axis=OUT),run_time=15,rate_func=linear)
        self.wait(1)        
        self.play(angle.animate.increment_value(-PI/2),Rotate(moon,-PI/2,axis=OUT),run_time=1)
        self.play(self.frame.animate.reorient(360,40,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3),run_time=1)
        self.play(sun.animate.become(Sun(big_glow_ratio=0.05).scale(0.01).move_to([RENDER_EARTH_RADIUS*6.5,0,0])),self.camera.light_source.animate.move_to([RENDER_EARTH_RADIUS*7.5,0,0]))   
        moon2= Moon().move_to([-RENDER_EARTH_RADIUS*6,0,0]).rotate(180*DEGREES,axis=OUT)
        moon2.move_to([-np.cos(220*DEGREES)*RENDER_EARTH_RADIUS*6,-np.sin(220*DEGREES)*RENDER_EARTH_RADIUS*6,0])
        self.play(FadeIn(moon2))
        self.wait()
        self.play(self.frame.animate.reorient(310,40,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3),run_time=1)
        self.wait(2)
        self.play(sun.animate.become(Sun(big_glow_ratio=0.05).scale(0.1).move_to([RENDER_EARTH_RADIUS*100,0,0])),self.camera.light_source.animate.move_to([RENDER_EARTH_RADIUS*100,0,0]))   
        self.wait(1)
        self.play(self.frame.animate.reorient(360,40,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3),run_time=1)
        self.embed()

class Parallax(Scene):
    def construct(self):
        # Starting with 20 arc-seconds
        tex1 = Tex(r"d_{(pc)}=\frac{1}{0,75''}").scale(2).shift(UP*2)
        tex2 = Tex(r"d_{(pc)}=\frac{1}{0,75''}",r"= 1,3 pc").scale(2).shift(UP*2)
        tex3 = Tex(r"41000000000000km").scale(2).shift(DOWN)
        tex4 = Tex(r"41000000000000km = 4,3 ly").scale(2).shift(DOWN)
        
        self.play(Write(tex1))
        self.wait(2)
        self.play(TransformMatchingTex(tex1, tex2),run_time=2)        
        self.wait(2)
        self.play(Write(tex3))
        self.wait(2)
        self.play(TransformMatchingTex(tex3, tex4),run_time=2)
        self.wait(2)

class EarthMoonSunSimilarTriangles_es(Scene):
    def construct(self):
        self.frame.reorient(-50, 80, 0, ORIGIN, 7)
        stars_data = extract_star_data()

        stars = Group()

        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            stars.add(star)
        self.add(stars)
        # Câmera
        d_moon = 12.0

        # Parâmetros fora de escala (somente para visualização)
        r_moon = 1.1
        d_sun  = 45.0
        r_sun  = 4.125

        # Cores
        col_earth = BLUE_E
        col_moon  = GREY_E
        col_sun   = YELLOW
        col_dist  = GREY_B
        col_rad   = GREY_B
        col_tri_m = TEAL_B
        col_tri_s = GOLD_E

        # Corpos celestes (com texturas originais do ManimGL/Saperientia)
        earth = Earth(clouds=False).scale(25).move_to(ORIGIN)
        moon  = Moon().scale(r_moon / 0.02).move_to(RIGHT * d_moon)
        sun   = Sun(big_glow_ratio=0.25).scale(r_sun / 7.0).move_to(RIGHT * d_sun)

        self.camera.light_source.move_to(sun.get_center())
        self.add(earth, moon, sun)
        self.wait(2)
        self.play(self.frame.animate.reorient(0, 0, 0, ORIGIN+RIGHT * d_moon, 22), run_time=2)
        # Rótulos
        # self.play(
        #     FadeIn(Tex("Terra").scale(2).next_to(earth, DOWN).set_color(col_earth)),
        #     FadeIn(Tex("Lua").scale(2).next_to(moon,  DOWN).set_color(col_moon)),
        #     FadeIn(Tex("Sol").scale(2).next_to(sun,   UP  ).set_color(col_sun))
        # )

        # Função: triângulo retângulo + marcador de 90° em “L”
        def right_triangle(index,origin, center, radius, tri_color):
            A = origin
            B = center
            C = B + UP * radius

            seg_AB = Line(A, B, color=col_dist, stroke_width=3).set_z_index(index)
            seg_BC = Line(B, C, color=col_rad,  stroke_width=3).set_z_index(index)
            seg_AC = Line(A, C, color=tri_color, stroke_width=4).set_z_index(index)

            # marcador de 90° feito manualmente
            size = 0.4
            corner = VGroup(
                Line(B + LEFT * 0.02+ UP    * size, B + LEFT * size+ UP    * size, color=tri_color, stroke_width=3),
                Line(B + UP    * 0.02+ LEFT * size, B + UP    * size+ LEFT * size, color=tri_color, stroke_width=3)
            )
            return VGroup(seg_AB, seg_BC, seg_AC, corner), A, B, C

        # Triângulos para Lua e Sol
        tri_moon, A_m, B_m, C_m = right_triangle(2,earth.get_center(), moon.get_center(), r_moon, col_tri_m)
        tri_sun,  A_s, B_s, C_s = right_triangle(5,earth.get_center(), sun.get_center(),  r_sun,  col_tri_s)

        self.play(ShowCreation(tri_moon),run_time=2)
        self.play(ShowCreation(tri_sun),run_time=2)

        # Rótulos das grandezas
        lbl_dm = Tex("d_{Luna}").set_color(col_tri_m).scale(1.5).next_to(Line(A_m, B_m).get_center(), DOWN)
        lbl_rm = Tex("r_{Luna}").set_color(col_tri_m).scale(1.5).next_to(Line(B_m, C_m).get_center(), RIGHT)
        lbl_ds = Tex("d_{Sol}").set_color(col_tri_s).scale(1.5).next_to(Line(A_s, B_s).get_center(), DOWN)
        lbl_rs = Tex("r_{Sol}").set_color(col_tri_s).scale(1.5).next_to(Line(B_s, C_s).get_center(), RIGHT)

        self.play(FadeIn(lbl_dm), FadeIn(lbl_rm),run_time=2)
        self.play(FadeIn(lbl_ds), FadeIn(lbl_rs),run_time=2)

        # Relação r/d
        ratio_m = Tex(r"\frac{r_{Luna}}{d_{Luna}}").set_color(col_tri_m).scale(0.9).scale(3)
        ratio_s = Tex(r"\frac{r_{Sol}}{d_{Sol}}").set_color(col_tri_s).scale(0.9).scale(3)
        approx  = Tex(r"\approx").scale(0.9).scale(3)
        VGroup(ratio_m, approx, ratio_s).arrange(RIGHT, buff=0.5).to_corner(UR).shift(UR*2)
        self.play(FadeIn(VGroup(ratio_m, approx, ratio_s)))

        # Guias verticais até o topo dos astros
        guide_moon = DashedLine(C_m, moon.get_center() + UP * r_moon, color=col_tri_m, stroke_width=2)
        guide_sun  = DashedLine(C_s, sun.get_center()  + UP * r_sun,  color=col_tri_s, stroke_width=2)
        self.play(ShowCreation(guide_moon))
        self.play(ShowCreation(guide_sun))

        self.play(self.frame.animate.reorient(0, 0, 0, ORIGIN+RIGHT * d_moon*2, 40), run_time=2)
        self.wait(3)
        
class EarthSun_VisualOnly(Scene):
    def construct(self):
        stars_data = extract_star_data()
        stars = Group(*[
            DotCloud(points=[np.array([x, y, z])],
                     color=color,
                     radius=size,
                     opacity=0.6)
            for x, y, z, size, color in stars_data
        ])
        self.add(stars)

        earth = Earth().scale(25)
        sun   = Sun()
        self.add(earth,sun)
        self.frame.reorient(90,90,0,ORIGIN,RENDER_SUN_RADIUS*3)
        earth.shift(RIGHT*15)
        self.wait(2)
        self.play(earth.animate.scale(1/25))
        self.play(self.frame.animate.reorient(80,90,0,earth.get_center(),RENDER_SUN_RADIUS*0.13),run_time=3)
        self.wait(0.5)
        self.frame.add_updater(lambda m:m.reorient(80,90,0,earth.get_center(),RENDER_SUN_RADIUS*0.2))
        self.play(earth.animate.move_to([RENDER_SUN_EARTH_DISTANCE,0,0]),run_time=3)



class Parallaxe(Scene):
    def construct(self):
        stars_data = extract_star_data()
        stars = Group(*[
            DotCloud(points=[np.array([x, y, z])],
                     color=color,
                     radius=size,
                     opacity=0.6)
            for x, y, z, size, color in stars_data
        ])
        self.add(stars)
        p=np.array([627751.2200000001, 6026669.100000001, -127126.83])
        self.frame.reorient(-10,+0.153788/PI*180+90,0,ORIGIN+[-9999,99999,0],10)
        self.play(self.frame.animate.reorient(-10,+0.153788/PI*180+90,0,ORIGIN+[99999,0,0],10),run_time=10)
        # self.frame.add_updater(lambda m:m.reorient(0,90,0,p,np.linalg.norm(p)))



class PenumbraAndUmbra(InteractiveScene):
    def construct(self):
        # Add earth and sun
        frame = self.frame
        frame.set_field_of_view(10 * DEG)
        light_source = self.camera.light_source

        sun = Sun(radius=2, big_glow_ratio=10, big_glow_factor=2)
        sun.move_to(7 * LEFT)

        earth = Earth(radius=0.7)
        earth.rotate(90 * DEG, LEFT).rotate(23*DEG, OUT)
        earth.move_to(2 * RIGHT)

        light_source.move_to(sun)
        self.add(sun, earth)

        # Add shadows
        umbra, penumbra, umbral_lines, penumbral_lines = shadow_group = self.get_umbral_lines(sun[0], earth)

        umbrum_word = Text("Umbra", font_size=24)
        umbrum_word.move_to(umbra)
        penumbrum_word = Text("Penumbra", font_size=24)
        penumbrum_word.move_to(umbrum_word)
        penumbrum_word.shift(0.55 * earth.get_height() * UP)

        self.add(penumbra, penumbral_lines, earth, penumbrum_word)
        self.play(
            FadeIn(penumbra),
            FadeIn(penumbral_lines),
            FadeIn(penumbrum_word),
        )
        self.add(*shadow_group, earth, umbrum_word, penumbrum_word)
        self.play(
            FadeIn(umbral_lines, lag_ratio=0),
            FadeIn(umbra, time_span=(0.5, 1)),
            FadeIn(umbrum_word, time_span=(0.5, 1)),
        )
        self.wait()

        # Pull away
        self.add(shadow_group, earth, umbrum_word, penumbrum_word)

        shift_vect = 400 * LEFT
        sun[2].set_z_index(-5)
        self.play(
            sun[:2].animate.shift(shift_vect),
            sun[2].animate.shift(5 * LEFT).scale(2),
            frame.animate.reorient(0, 0, 0, (4.21, 0.56, 0.0), 13.04),
            umbrum_word.animate.shift(2 * RIGHT),
            penumbrum_word.animate.scale(0.5, about_edge=DOWN).shift(2 * RIGHT),
            UpdateFromFunc(shadow_group, lambda m: m.become(self.get_umbral_lines(sun[0], earth))),
            run_time=8,
        )
        moon = Moon(radius=0.7/4).move_to(umbra.get_center()*0.08+DOWN*3*0.7/4).set_opacity(0.7)
        moon2 = Moon(radius=0.7/4).move_to(umbra.get_center()*0.08+DOWN*1*0.7/4).set_opacity(0.7)
        moon3 = Moon(radius=0.7/4).move_to(umbra.get_center()*0.08+UP*1*0.7/4).set_opacity(0.7)
        moon4 = Moon(radius=0.7/4).move_to(umbra.get_center()*0.08+UP*3*0.7/4).set_opacity(0.7)
        self.play(FadeIn(moon))
        self.play(FadeIn(moon2))
        self.play(FadeIn(moon3))
        self.play(FadeIn(moon4))


    def get_umbral_lines(self, circle1, circle2):
        r1 = circle1.get_width() / 2
        r2 = circle2.get_width() / 2
        c1 = circle1.get_center()
        c2 = circle2.get_center()

        vect = c2 - c1
        dist = get_norm(vect)

        cs1 = c1 + vect / (1.0 - r2 / r1)
        cs2 = c1 + vect / (1.0 + r2 / r1)

        angle = math.asin(r1 / get_norm(cs1 - c1))

        v1 = rotate_vector(UP, -angle)
        v2 = rotate_vector(DOWN, angle)
        v3 = rotate_vector(UP, angle)
        v4 = rotate_vector(DOWN, -angle)

        umbral_lines = VGroup(
            Line(c1 + r1 * v1, cs1),
            Line(c1 + r1 * v2, cs1),
        )
        umbral_lines.set_stroke(WHITE, 0.5)
        umbra = Polygon(cs1, c2 + r2 * v1, c2, c2 + r2 * v2)
        umbra.set_fill(BLACK, 1).set_stroke(width=0)

        penumbral_lines = VGroup(
            Line(cs2, c2 + r2 * v3).scale(10),
            Line(cs2, c2 + r2 * v4).scale(10),
        )
        penumbral_lines.set_stroke(WHITE, 0.5)
        penumbra = Polygon(
            c2 + r2 * v3, penumbral_lines[0].get_end(),
            penumbral_lines[1].get_end(), c2 + r2 * v4,
        )
        penumbra.set_fill(BLACK, 0.5)
        penumbra.set_stroke(width=0)

        return VGroup(umbra, penumbra, umbral_lines, penumbral_lines)
    
class NearbyStarsPT(InteractiveScene):
    def construct(self):
        # Add sun and earth
        orbit_radius = 3.5
        conversion_factor = orbit_radius / SUN_EARTH_DISTANCE

        sun = Sun(radius=conversion_factor * SUN_RADIUS, big_glow_ratio=20)
        sun.center()
        orbit = Circle(radius=orbit_radius)
        orbit.set_stroke(BLUE, (0, 4))
        earth_glow = GlowDot(color=BLUE)
        earth_glow.f_always.move_to(orbit.get_start)


        self.add( sun, orbit, earth_glow)

        # Show the astronomical unit
        dist_line = Line()
        dist_line.set_stroke(WHITE, 1)
        dist_line.f_always.put_start_and_end_on(sun.get_center, orbit.get_start)

        dist_label = Text("Unidade\nAstronômica", font_size=36)
        dist_label.f_always.move_to(
            lambda: dist_line.get_center() + 0.5 * normalize(rotate_vector(dist_line.get_vector(), 90 * DEG))
        )

        self.play(
            FadeIn(dist_line, time_span=(0, 1)),
            FadeIn(dist_label, time_span=(0, 1)),
            Rotate(orbit, TAU, about_point=ORIGIN, rate_func=linear, run_time=10),
        )
        self.wait()

        # Transition to initials
        dist_label.clear_updaters()
        au_label = Text("U.A.", font_size=36)

        def update_au_label(label):
            point = dist_line.get_center()
            direction = normalize(rotate_vector(point, 90 * DEG))
            step = 0.65 * interpolate(label.get_width(), label.get_height(), abs(direction[1]))
            label.move_to(point + step * direction)

        au_label.add_updater(update_au_label)

        self.play(LaggedStart(
            *(
                ReplacementTransform(dist_label[t2][0], au_label[t1][i])
                for t1, t2, i in zip("U.A.", ["U", "nidade", "A", "stronômica"], [0, 0, 0, 1])
            ),
            lag_ratio=0.2
        ))
        self.add(au_label)

        # Position to the side
        frame = self.frame
        self.play(
            Rotate(orbit, 90 * DEG),
            frame.animate.reorient(0, 0, 0, 7 * RIGHT, 14),
            run_time=2
        )

        # Zoom into and out of earth real quick
        frame.save_state()
        earth = Earth(radius=orbit_radius * (EARTH_RADIUS / SUN_EARTH_DISTANCE))
        earth.move_to(earth_glow)
        earth.rotate(23*DEG, RIGHT)
        frame.move_to(earth)
        frame.set_height(2 * earth.get_height())
        frame.reorient(-74, 79, 0)
        self.camera.light_source.move_to(sun)

        self.remove(earth_glow, orbit, dist_line)
        self.add(earth)
        self.wait()
        srf = squish_rate_func(smooth, 0.7, 1)
        self.play(
            UpdateFromAlphaFunc(frame, lambda m, a: m.reorient(
                *interpolate(np.array([-74, 79, 0]), np.zeros(3), a),
                interpolate(earth.get_center(), 7 * RIGHT, srf(a)),
                np.exp(interpolate(np.log(2 * earth.get_height()), np.log(14), smooth(a))),
            ), run_time=5),
            FadeIn(earth_glow, time_span=(2.5, 4.5)),
            FadeIn(orbit, time_span=(1, 4)),
            FadeIn(dist_line, time_span=(1, 4)),
            FadeIn(au_label, time_span=(4, 5)),
            FadeOut(earth),
            run_time=5,
        )

        # Show observations
        star = Group(
            ImageMobject('https://images.vexels.com/media/users/3/254382/isolated/preview/8efce08800d999b79c2f73b94c75fd03-estrela-amarela-plana.png').set_height(0.8).center(),
            GlowDot(color=WHITE).center()
        )
        star[1].add_updater(lambda m: m.set_width(0.4 * ((1 + math.sin(1.5 * self.time)))))
        star.move_to(50 * RIGHT)
        obs_points = Group(
            TrueDot(point, radius=0.1).set_color(GREEN).make_3d()
            for point in [orbit.get_top(), orbit.get_bottom()]
        )
        obs_lines = VGroup(
            self.get_obs_line(obs_point, star)
            for obs_point in obs_points
        )
        obs_lines.set_stroke(WHITE, 2)
        for line, point in zip(obs_lines, obs_points):
            line.start_point = point
            line.star = star
            line.add_updater(lambda m: m.put_start_and_end_on(m.start_point.get_center(), m.star.get_center()))

        obs_labels = VGroup(Text(f"Observação {n}") for n in [1, 2])
        for label, point, vect in zip(obs_labels, obs_points, [UP, DOWN]):
            label.next_to(point, vect, MED_SMALL_BUFF)

        self.add(star)

        self.play(
            ShowCreation(obs_lines[0], suspend_mobject_updating=True),
            FadeIn(obs_labels[0], 0.25 * UP),
            FadeIn(obs_points[0]),
        )
        self.wait()
        self.play(Rotate(orbit, PI), run_time=2)
        self.play(
            ShowCreation(obs_lines[1], suspend_mobject_updating=True),
            FadeIn(obs_labels[1], DOWN),
            FadeIn(obs_points[1]),
        )
        self.wait()

        # Show the angle vary during the orbit
        self.play(
            star.animate.move_to(15 * RIGHT),
            run_time=2
        )
        self.wait()

        obs_lines.suspend_updating()
        sample_obs_line = self.get_obs_line(earth_glow, star)
        self.play(
            FadeIn(sample_obs_line),
            obs_lines.animate.set_stroke(opacity=0.1)
        )
        self.play(Rotate(orbit, PI, run_time=10))
        self.wait()
        self.play(
            FadeOut(sample_obs_line),
            obs_lines.animate.set_stroke(opacity=1),
        )

        # # Pull it far away, then back
        # curr_center = star.get_center()
        curr_angle = obs_lines[1].get_angle() - obs_lines[0].get_angle()
        # orbit_radius / math.tan(curr_angle / 2)

        # obs_lines.resume_updating()
        # self.play(
        #     UpdateFromAlphaFunc(star, lambda m, a: m.move_to(
        #         RIGHT * orbit_radius / math.tan(interpolate(curr_angle, 1e-5, there_and_back_with_pause(a)) / 2)
        #     )),
        #     run_time=6,
        # )

        # Label the distance and angle
        line_to_star = Line(sun.get_center(), star.get_center())
        line_to_star.set_stroke(RED, 3)
        dist_label = Tex("D", font_size=60)
        dist_label.next_to(line_to_star, UP, buff=2 * SMALL_BUFF)
        dist_label.match_color(line_to_star)

        arc = Arc(PI, -curr_angle / 2, arc_center=star.get_center(), radius=3)
        arc_label = Tex(R"\theta / 2", font_size=60)
        arc_label.next_to(arc, LEFT, buff=SMALL_BUFF)

        self.play(
            ShowCreation(line_to_star),
            obs_lines.animate.set_stroke(width=1),
            FadeIn(dist_label, RIGHT),
        )
        self.wait()
        self.play(
            ShowCreation(arc),
            Write(arc_label),
        )
        self.play(FlashAround(arc_label, run_time=2))
        self.wait()
        self.play(
            Transform(obs_lines[0].copy().clear_updaters(), obs_lines[1].copy(), remover=True),
            run_time=2
        )
        self.wait()

        #

    def get_obs_line(self, obj1, obj2, dash_length=0.1, stroke_color=WHITE, stroke_width=2):
        # line = DashedLine(obj1.get_center(), obj2.get_center())
        line = Line(obj1.get_center(), obj2.get_center())
        line.set_stroke(stroke_color, stroke_width)
        line.f_always.put_start_and_end_on(obj1.get_center, obj2.get_center)
        return line


class NearbyStarsEs(InteractiveScene):
    def construct(self):
        # Add sun and earth
        orbit_radius = 3.5
        conversion_factor = orbit_radius / SUN_EARTH_DISTANCE

        sun = Sun(radius=conversion_factor * SUN_RADIUS, big_glow_ratio=20)
        sun.center()
        orbit = Circle(radius=orbit_radius)
        orbit.set_stroke(BLUE, (0, 4))
        earth_glow = GlowDot(color=BLUE)
        earth_glow.f_always.move_to(orbit.get_start)


        self.add( sun, orbit, earth_glow)

        # Show the astronomical unit
        dist_line = Line()
        dist_line.set_stroke(WHITE, 1)
        dist_line.f_always.put_start_and_end_on(sun.get_center, orbit.get_start)

        dist_label = Text("Unidad\nAstronómica", font_size=36)
        dist_label.f_always.move_to(
            lambda: dist_line.get_center() + 0.5 * normalize(rotate_vector(dist_line.get_vector(), 90 * DEG))
        )

        self.play(
            FadeIn(dist_line, time_span=(0, 1)),
            FadeIn(dist_label, time_span=(0, 1)),
            Rotate(orbit, TAU, about_point=ORIGIN, rate_func=linear, run_time=10),
        )
        self.wait()

        # Transition to initials
        dist_label.clear_updaters()
        au_label = Text("U.A.", font_size=36)

        def update_au_label(label):
            point = dist_line.get_center()
            direction = normalize(rotate_vector(point, 90 * DEG))
            step = 0.65 * interpolate(label.get_width(), label.get_height(), abs(direction[1]))
            label.move_to(point + step * direction)

        au_label.add_updater(update_au_label)

        self.play(LaggedStart(
            *(
                ReplacementTransform(dist_label[t2][0], au_label[t1][i])
                for t1, t2, i in zip("U.A.", ["U", "nidad", "A", "stronómica"], [0, 0, 0, 1])
            ),
            lag_ratio=0.2
        ))
        self.add(au_label)

        # Position to the side
        frame = self.frame
        self.play(
            Rotate(orbit, 90 * DEG),
            frame.animate.reorient(0, 0, 0, 7 * RIGHT, 14),
            run_time=2
        )

        # Zoom into and out of earth real quick
        frame.save_state()
        earth = Earth(radius=orbit_radius * (EARTH_RADIUS / SUN_EARTH_DISTANCE))
        earth.move_to(earth_glow)
        earth.rotate(23*DEG, RIGHT)
        frame.move_to(earth)
        frame.set_height(2 * earth.get_height())
        frame.reorient(-74, 79, 0)
        self.camera.light_source.move_to(sun)

        self.remove(earth_glow, orbit, dist_line)
        self.add(earth)
        self.wait()
        srf = squish_rate_func(smooth, 0.7, 1)
        self.play(
            UpdateFromAlphaFunc(frame, lambda m, a: m.reorient(
                *interpolate(np.array([-74, 79, 0]), np.zeros(3), a),
                interpolate(earth.get_center(), 7 * RIGHT, srf(a)),
                np.exp(interpolate(np.log(2 * earth.get_height()), np.log(14), smooth(a))),
            ), run_time=5),
            FadeIn(earth_glow, time_span=(2.5, 4.5)),
            FadeIn(orbit, time_span=(1, 4)),
            FadeIn(dist_line, time_span=(1, 4)),
            FadeIn(au_label, time_span=(4, 5)),
            FadeOut(earth),
            run_time=5,
        )

        # Show observations
        star = Group(
            ImageMobject('https://images.vexels.com/media/users/3/254382/isolated/preview/8efce08800d999b79c2f73b94c75fd03-estrela-amarela-plana.png').set_height(0.8).center(),
            GlowDot(color=WHITE).center()
        )
        star[1].add_updater(lambda m: m.set_width(0.4 * ((1 + math.sin(1.5 * self.time)))))
        star.move_to(50 * RIGHT)
        obs_points = Group(
            TrueDot(point, radius=0.1).set_color(GREEN).make_3d()
            for point in [orbit.get_top(), orbit.get_bottom()]
        )
        obs_lines = VGroup(
            self.get_obs_line(obs_point, star)
            for obs_point in obs_points
        )
        obs_lines.set_stroke(WHITE, 2)
        for line, point in zip(obs_lines, obs_points):
            line.start_point = point
            line.star = star
            line.add_updater(lambda m: m.put_start_and_end_on(m.start_point.get_center(), m.star.get_center()))

        obs_labels = VGroup(Text(f"Observación {n}") for n in [1, 2])
        for label, point, vect in zip(obs_labels, obs_points, [UP, DOWN]):
            label.next_to(point, vect, MED_SMALL_BUFF)

        self.add(star)

        self.play(
            ShowCreation(obs_lines[0], suspend_mobject_updating=True),
            FadeIn(obs_labels[0], 0.25 * UP),
            FadeIn(obs_points[0]),
        )
        self.wait()
        self.play(Rotate(orbit, PI), run_time=2)
        self.play(
            ShowCreation(obs_lines[1], suspend_mobject_updating=True),
            FadeIn(obs_labels[1], DOWN),
            FadeIn(obs_points[1]),
        )
        self.wait()

        # Show the angle vary during the orbit
        self.play(
            star.animate.move_to(15 * RIGHT),
            run_time=2
        )
        self.wait()

        obs_lines.suspend_updating()
        sample_obs_line = self.get_obs_line(earth_glow, star)
        self.play(
            FadeIn(sample_obs_line),
            obs_lines.animate.set_stroke(opacity=0.1)
        )
        self.play(Rotate(orbit, PI, run_time=10))
        self.wait()
        self.play(
            FadeOut(sample_obs_line),
            obs_lines.animate.set_stroke(opacity=1),
        )

        # # Pull it far away, then back
        # curr_center = star.get_center()
        curr_angle = obs_lines[1].get_angle() - obs_lines[0].get_angle()
        # orbit_radius / math.tan(curr_angle / 2)

        # obs_lines.resume_updating()
        # self.play(
        #     UpdateFromAlphaFunc(star, lambda m, a: m.move_to(
        #         RIGHT * orbit_radius / math.tan(interpolate(curr_angle, 1e-5, there_and_back_with_pause(a)) / 2)
        #     )),
        #     run_time=6,
        # )

        # Label the distance and angle
        line_to_star = Line(sun.get_center(), star.get_center())
        line_to_star.set_stroke(RED, 3)
        dist_label = Tex("D", font_size=60)
        dist_label.next_to(line_to_star, UP, buff=2 * SMALL_BUFF)
        dist_label.match_color(line_to_star)

        arc = Arc(PI, -curr_angle / 2, arc_center=star.get_center(), radius=3)
        arc_label = Tex(R"\theta / 2", font_size=60)
        arc_label.next_to(arc, LEFT, buff=SMALL_BUFF)

        self.play(
            ShowCreation(line_to_star),
            obs_lines.animate.set_stroke(width=1),
            FadeIn(dist_label, RIGHT),
        )
        self.wait()
        self.play(
            ShowCreation(arc),
            Write(arc_label),
        )
        self.play(FlashAround(arc_label, run_time=2))
        self.wait()
        self.play(
            Transform(obs_lines[0].copy().clear_updaters(), obs_lines[1].copy(), remover=True),
            run_time=2
        )
        self.wait()

        #

    def get_obs_line(self, obj1, obj2, dash_length=0.1, stroke_color=WHITE, stroke_width=2):
        # line = DashedLine(obj1.get_center(), obj2.get_center())
        line = Line(obj1.get_center(), obj2.get_center())
        line.set_stroke(stroke_color, stroke_width)
        line.f_always.put_start_and_end_on(obj1.get_center, obj2.get_center)
        return line



def get_sphere_mesh(radius=1.0):
    sphere = Sphere(radius=radius)
    mesh = SurfaceMesh(sphere)
    mesh.set_stroke(WHITE, 0.5, 0.5)
    return mesh

EARTH_TILT_ANGLE =0*DEG

class SizeOfEarthRenewed(InteractiveScene):
    radius = 3.0

    def construct(self):
        # Setup
        self.frame.reorient(center=RIGHT*3,height=10)
        self.set_floor_plane("xz")
        frame = self.frame
        frame.set_field_of_view(30 * DEG)

        light = self.camera.light_source
        light.move_to(20 * RIGHT)

        # Add earth
        sphere = Sphere(radius=self.radius)
        earth = Earth(radius=self.radius).rotate(-260*DEG,OUT).rotate(-90*DEG,RIGHT).set_z_index(2)




        slice_tracker = ValueTracker(self.radius + SMALL_BUFF)
        earth.add_updater(lambda m: m.set_clip_plane(OUT, slice_tracker.get_value()))


        self.add(earth)

        # Unflatten earth
        earth.save_state()
        earth.data["d_normal_point"] = earth.get_points() + 1e-3 * RIGHT
        earth.note_changed_data()



        # Add rays from the sun
        sun = GlowDot(5 * RIGHT, radius=1)
        n_rays = 25
        rays = Line(LEFT, RIGHT).replicate(n_rays)
        rays.set_stroke(YELLOW, 1.5)

        def update_rays(rays):
            ys = np.linspace(earth.get_y(UP), earth.get_y(DOWN), len(rays))
            for ray, y in zip(rays, ys):
                ray.put_start_and_end_on(
                    sun.get_center(),
                    [math.sqrt(abs(self.radius**2 - y**2)), y, 0],
                )

        rays.add_updater(update_rays)
        
        n_rays2 = 50
        rays2 = Line(LEFT, RIGHT).replicate(n_rays2)
        rays2.set_stroke(YELLOW, 1.5)

        def update_rays2(rays):
            ys = np.linspace(0, 2 * np.pi, n_rays2, endpoint=False)
            for ray, y in zip(rays, ys):
                ray.put_start_and_end_on(
                    sun.get_center(),
                    [math.cos(y)*2000, math.sin(y)*2000, 0]+sun.get_center(),
                )

        rays2.add_updater(update_rays2,False)
        self.play(FadeIn(sun), run_time=2)
        self.wait(3)
        self.play(FadeIn(rays2))
        self.wait(3)
        self.play(FadeIn(rays),rays2.animate.set_opacity(0.3))
        self.wait(1)
        
        
        
        


        self.wait(5)
        self.play(
            sun.animate.move_to(1000 * RIGHT),
            run_time=10
        )

class Raios(Scene):
    def construct(self):
        self.frame.reorient(90,70,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3)
        esfera_celeste=EsferaCeleste()
        esfera_celeste.add_updater(lambda m: self.add(esfera_celeste))
        superficie = SuperficieObservador()
        self.add(superficie)
        self.add(esfera_celeste)
        sun = Sun(radius=RENDER_CELESTIAL_SPHERE_RADIUS*0.08).move_to(P(48,180).get_center())
        self.frame.add_ambient_rotation(4*DEG)
        self.play(FadeIn(sun))
        o = Line3D(ORIGIN,ORIGIN+[0,0,CELESTIAL_SPHERE_RADIUS*0.2],color=PURPLE,width=0.003)
        self.play(FadeIn(o))
        
        
        n_rays2 = 6
        rays = [Line3D(LEFT, RIGHT,width=0.03,color=YELLOW) for i in range(0,n_rays2)]
        group = [rays for i in range(0,n_rays2)]
        raios=Group()
        factor = 0.06
        for n, ray in enumerate(group):
            for i,r in enumerate(ray):
                r=Line3D(LEFT, RIGHT,width=0.01,color=YELLOW).set_opacity(0.7)
                r.put_start_and_end_on(
                    (n-(n_rays2-1)/2)*RIGHT*factor+(i-(n_rays2-1)/2)*DOWN*factor,
                    (n-(n_rays2-1)/2)*RIGHT*factor+0.2*OUT+(i-(n_rays2-1)/2)*DOWN*factor+3*DOWN*factor,
                )
                raios.add(r)
        self.play(FadeIn(raios))
        self.wait(3)
        linha_astro = MarcadorAltura(P(48,180),arco=False,espessura_linha=5,cor_linha=RED)
        self.play(FadeIn(linha_astro))
        self.wait(5)
        self.frame.add_ambient_rotation(-4*DEG)
        self.play(self.frame.animate.reorient(0,90-48,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3))
        self.wait(2)
        self.play(self.frame.animate.reorient(360-160,70,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3))
        
        self.wait(2)
        self.frame.add_ambient_rotation(-4*DEG)
        o2 = Line3D(ORIGIN,ORIGIN+[0,0,CELESTIAL_SPHERE_RADIUS*0.2],color=BLUE,width=0.003).shift(LEFT*0.12)
        self.play(FadeIn(o2))
        self.play(raios.animate.shift(RIGHT*0.12),o.animate.shift(RIGHT*0.12),o2.animate.shift(RIGHT*0.12))
        
        esfera_celeste2=EsferaCeleste().shift(RIGHT*0.12)
        self.play(FadeIn(esfera_celeste2))
        esfera_celeste2.add_updater(lambda m: self.add(esfera_celeste2))
        superficie2 = SuperficieObservador().shift(RIGHT*0.12)
        self.play(FadeIn(superficie2))
        sun2 = Sun(radius=RENDER_CELESTIAL_SPHERE_RADIUS*0.08).move_to(P(48,180).get_center()).shift(RIGHT*0.12)
        self.play(FadeIn(sun2))
        linha_astro2 = MarcadorAltura(P(48,180),arco=False,espessura_linha=5,cor_linha=RED).shift(RIGHT*0.12)
        self.play(FadeIn(linha_astro2))
        self.wait(4)     


class Qst2(Scene):
    def construct(self):

        # Define the values
        parsec_to_au = 206265
        au_to_meters = 150_000_000_000  # 1.5 x 10^11 meters

        # Calculate 1 parsec in meters
        parsec_in_meters = parsec_to_au * au_to_meters

        # Create Manim objects for display
        text1 = Tex(r"1 \mathrm{parsec} = 206265 \mathrm{ U.A.}").scale(1.5).shift(UP)
        text2 = Tex(r"1 \mathrm{U.A.} = 1.5 \times 10^{11} \mathrm{ m}").scale(1.5).next_to(text1, DOWN) # This line is already correct
        
        # Display the calculation
        calc_text = Tex(
            r"1 \mathrm{ parsec} = 206265 \times (1.5 \times 10^{11}) \mathrm{ m}"
        ).scale(1.5).next_to(text2, 2*DOWN)

        # Display the result in scientific notation
        res = Tex(
            r"1 \mathrm{ parsec} = 3,09 \times 10^{15} \mathrm{ m}"
        ).scale(1.5).next_to(calc_text, DOWN)


        self.play(Write(text1))
        self.wait(1)
        self.play(Write(text2))
        self.wait(1)
        self.play(Write(calc_text))
        self.wait(2)
        self.play(Write(res))
        self.wait(3)

class LightYearCalculation(Scene):
    def construct(self):
        # Python values for calculation
        c = 3e8
        seconds_in_year = 365.25 * 24 * 60 * 60
        light_year_in_meters = c * seconds_in_year

        # 1. Display speed of light
        light_speed_tex = Tex(r"c = 3 \times 10^8 \mathrm{ m/s}").scale(1.5).to_edge(UP)
        
        self.play(Write(light_speed_tex))
        self.wait(1)

        # 2. Display calculation for seconds in a year
        seconds_calc_tex = Tex(
            r"1 \mathrm{ ano} = 365,25 \mathrm{ dias} \times 24 \mathrm{ h} \times 60 \mathrm{ min} \times 60 \mathrm{ s}"
        ).scale(1.2).next_to(light_speed_tex, DOWN, buff=LARGE_BUFF)
        
        self.play(Write(seconds_calc_tex))
        self.wait(2)

        # 3. Display result for seconds in a year
        seconds_res_tex = Tex(
            r"1 \mathrm{ ano} = 31 557 600 \mathrm{ s}"
        ).scale(1.2).next_to(seconds_calc_tex, DOWN, buff=MED_LARGE_BUFF)
        
        self.play(TransformFromCopy(seconds_calc_tex, seconds_res_tex))
        self.wait(2)

        # 4. Display final calculation for light-year
        light_year_calc_tex = Tex(
            r"1 \mathrm{ ly} = (3 \times 10^8 \mathrm{ m/s}) \times (31557600 \mathrm{ s})"
        ).scale(1.2).next_to(seconds_res_tex, DOWN, buff=LARGE_BUFF)
        
        self.play(Write(light_year_calc_tex))
        self.wait(2)

        # 5. Display final result
        final_res_tex = Tex(
            r"1 \mathrm{ ly} \approx 9,46 \times 10^{15} \mathrm{ m}"
        ).scale(1.5).next_to(light_year_calc_tex, DOWN, buff=MED_LARGE_BUFF)
        
        self.play(TransformFromCopy(light_year_calc_tex, final_res_tex))
        self.wait(3)

        # 6. Fade out old calculations and reposition the light-year value
        self.play(
            FadeOut(light_speed_tex),
            FadeOut(seconds_calc_tex),
            FadeOut(seconds_res_tex),
            FadeOut(light_year_calc_tex),
        )
        self.wait(0.5)
        self.play(final_res_tex.animate.to_edge(UP))
        self.wait(1)

        # 7. Display the value of 1 parsec in meters
        parsec_in_meters_tex = Tex(r"1 \mathrm{ parsec} \approx 3,09 \times 10^{16} \mathrm{ m}").scale(1.5)
        parsec_in_meters_tex.next_to(final_res_tex, DOWN, buff=MED_LARGE_BUFF)
        
        self.play(Write(parsec_in_meters_tex))
        self.wait(2)

        # 8. Display the division to find parsecs in light-years
        division_calc_tex = Tex(
            r"1 \mathrm{ parsec} = \frac{3,09 \times 10^{16} \mathrm{ m}}{9,46 \times 10^{15} \mathrm{ m}}"
        ).scale(1.5).next_to(parsec_in_meters_tex, DOWN, buff=LARGE_BUFF)
        
        self.play(Write(division_calc_tex))
        self.wait(2)

        # 9. Display the final result of the conversion
        final_conversion_res_tex = Tex(
            r"1 \mathrm{ parsec} \approx 3,26 \mathrm{ ly}"
        ).scale(1.5).next_to(division_calc_tex, DOWN, buff=MED_LARGE_BUFF)
        
        self.play(TransformFromCopy(division_calc_tex, final_conversion_res_tex))
        self.wait(3)

class ParallaxInverseRelationship(Scene):
    def construct(self):
        # 1. Display the formula
        formula = Tex(r"d_{(\mathrm{pc})} = \frac{1}{p_{('')}}").scale(2)
        self.play(Write(formula))
        self.wait(2)
        self.play(formula.animate.to_edge(UP))
        self.wait(1)

        # 2. Configura os trackers e a exibição da equação
        p_tracker = ValueTracker(1.0)

        # Cria as partes da equação
        d_label = Tex(r"d_{(\mathrm{pc})}=").scale(1.5)
        one = Tex("1").scale(1.5)
        p_value = DecimalNumber(p_tracker.get_value(), num_decimal_places=3).scale(1.5)
        div_line = Line(LEFT, RIGHT)

        fraction = VGroup(one, div_line, p_value)

        equals_sign = Tex("=").scale(1.5)
        d_value = DecimalNumber(1.0, num_decimal_places=0).scale(1.5)
        
        d_units = Tex(r"\,\mathrm{pc}").scale(1.5)

        # Agrupa as partes da equação
        full_equation = VGroup(d_label, fraction, equals_sign, d_value, d_units)
        full_equation.next_to(formula, DOWN, buff=LARGE_BUFF)

        # Adiciona updaters para manter tudo alinhado
        p_value.add_updater(lambda m: m.set_value(p_tracker.get_value()))
        d_value.add_updater(lambda m: m.set_value(
            1 / p_tracker.get_value() 
        ))

        def layout_updater(m):
            new_width = max(one.get_width(), p_value.get_width()) + 0.1
            div_line.set_width(new_width)
            fraction.arrange(DOWN, buff=0.35)
            m.arrange(RIGHT, buff=MED_SMALL_BUFF)
            m.next_to(formula, DOWN, buff=LARGE_BUFF)

        full_equation.add_updater(layout_updater)

        self.play(
            Write(full_equation)
        )
        self.wait(2)

        # 3. Anima a diminuição de p, mostrando o aumento de d
        self.play(p_tracker.animate.set_value(0.1), run_time=3, rate_func=smooth)
        self.wait(1)
        self.play(p_tracker.animate.set_value(0.01), run_time=3, rate_func=smooth)
        self.wait(2)

        # 4. Mostra o limite
        self.play(p_tracker.animate.set_value(0.001), run_time=2, rate_func=smooth)
        self.wait(1)
        p_value2=DecimalNumber(p_tracker.get_value(), num_decimal_places=6).scale(1.5)
        self.play(Transform(p_value,p_value2))
        p_value2.add_updater(lambda m: m.set_value(p_tracker.get_value()))
        self.wait(1)
        self.play(p_tracker.animate.set_value(0.000000001), run_time=2, rate_func=smooth)
        self.wait(1)
        # Limpa os updaters antes de transformar
        full_equation.clear_updaters()
        p_value.clear_updaters()
        d_value.clear_updaters()

        final_text = Tex(r"d \to \infty").scale(2).move_to(full_equation)

        self.play(
            ReplacementTransform(full_equation, final_text)
        )
        self.wait(3)

class ParallaxCalculationExample(Scene):
    def construct(self):
        # 1. Display the formula
        formula = Tex(r"d_{(\mathrm{pc})} = \frac{1}{p_{('')}}").scale(1.5).to_edge(UP)

        # 2. Display the given parallax value in milliarcseconds
        p_mas = Tex(r"p = 5,4 \, \mathrm{mas}").scale(1.5).next_to(formula, DOWN, buff=LARGE_BUFF)
        
        # 3. Parallax in arcseconds
        p_arcsec = Tex(r"p = 0,0054''").scale(1.5).move_to(p_mas)

        # 4. Substitute into the formula
        calculation = Tex(r"d = \frac{1}{0,0054}").scale(1.5).next_to(p_arcsec, DOWN, buff=LARGE_BUFF)

        # 5. Show the final result
        result = Tex(r"d \approx 185,18 \, \mathrm{pc}").scale(1.5).next_to(calculation, DOWN)

        # Animation
        self.play(Write(formula))
        self.wait(1)
        
        self.play(Write(p_mas))
        self.wait(1)
        
        self.play(ReplacementTransform(p_mas, p_arcsec))
        self.wait(2)
        
        self.play(TransformFromCopy(VGroup(formula, p_arcsec), calculation))
        self.wait(2)
        
        self.play(TransformFromCopy(calculation, result))
        self.wait(3)

class ParallaxInverseProportionality(Scene):
    def construct(self):
        # Define a escala para o texto
        scale_factor = 1.8

        # 1. Mostra d = 1/p
        formula1 = Tex(r"d = \frac{1}{p}").scale(scale_factor)
        self.play(Write(formula1))
        self.wait(1)

        # 2. Transforma para p = 1/d
        formula2 = Tex(r"p = \frac{1}{d}").scale(scale_factor)
        self.play(TransformMatchingTex(formula1, formula2, path_arc=-PI/2))
        self.wait(2)

        # 3. Mostra o que acontece quando d -> 2d
        formula2_copy = formula2.copy()
        self.play(formula2.animate.shift(UP*2))
        step1 = Tex(r"p' = \frac{1}{2d}").scale(scale_factor)
        self.play(TransformMatchingTex(formula2_copy, step1, key_map={"p": "p'", "d": "2d"}))
        self.wait(2)

        # 4. Reorganiza para mostrar a relação
        step2 = Tex(r"p' = \frac{1}{2} \times \frac{1}{d}").scale(scale_factor)
        self.play(TransformMatchingTex(step1, step2))
        self.wait(2)

        # 5. Como p = 1/d
        step3 = Tex(r"p' = \frac{1}{2} \times p").scale(scale_factor)
        self.play(TransformMatchingTex(step2, step3, key_map={r"\frac{1}{d}": " p"}))
        self.wait(2)

        # 6. Forma final
        step4 = Tex(r"p' = \frac{p}{2}").scale(scale_factor)
        self.play(TransformMatchingTex(step3, step4,key_map={r"1": " p"}))
        self.wait(2)

        # 7. Conclusão final
        conclusion_group = VGroup(
            Tex(r"d \rightarrow 2d").scale(scale_factor),
            Tex(r"\Rightarrow").scale(scale_factor),
            Tex(r"p \rightarrow \frac{p}{2}").scale(scale_factor)
        ).arrange(RIGHT, buff=MED_LARGE_BUFF)

        self.play(
            FadeOut(formula2),
            ReplacementTransform(step4, conclusion_group)
        )
        self.wait(3)

class NoInfinito(InteractiveScene):
    def construct(self):
        # Add sun and earth
        orbit_radius = 3.5
        conversion_factor = orbit_radius / SUN_EARTH_DISTANCE

        sun = Sun(radius=conversion_factor * SUN_RADIUS, big_glow_ratio=20)
        sun.center()
        orbit = Circle(radius=orbit_radius)
        orbit.set_stroke(BLUE, (0, 4))
        earth_glow = GlowDot(color=BLUE)
        earth_glow.f_always.move_to(orbit.get_start)


        self.add( sun, orbit, earth_glow)


        # Position to the side
        frame = self.frame
        self.play(
            Rotate(orbit, 90 * DEG),
            run_time=2
        )

        # Zoom into and out of earth real quick
        frame.save_state()
        earth = Earth(radius=orbit_radius * (EARTH_RADIUS / SUN_EARTH_DISTANCE))
        earth.move_to(earth_glow)
        earth.rotate(23*DEG, RIGHT)
        frame.move_to(earth)
        frame.set_height(2 * earth.get_height())
        frame.reorient(-74, 79, 0)
        self.camera.light_source.move_to(sun)

        self.remove(earth_glow, orbit)
        self.add(earth)
        self.wait()
        srf = squish_rate_func(smooth, 0.7, 1)
        self.play(
            UpdateFromAlphaFunc(frame, lambda m, a: m.reorient(
                *interpolate(np.array([-74, 79, 0]), np.zeros(3), a),
                interpolate(earth.get_center(), 7 * RIGHT, srf(a)),
                np.exp(interpolate(np.log(2 * earth.get_height()), np.log(14), smooth(a))),
            ), run_time=5),
            FadeIn(earth_glow, time_span=(2.5, 4.5)),
            FadeIn(orbit, time_span=(1, 4)),
            FadeOut(earth),
            run_time=5,
        )

        # Show observations
        star = Group(
            ImageMobject('https://images.vexels.com/media/users/3/254382/isolated/preview/8efce08800d999b79c2f73b94c75fd03-estrela-amarela-plana.png').set_height(0.8).center(),
            GlowDot(color=WHITE).center()
        )
        star[1].add_updater(lambda m: m.set_width(0.4 * ((1 + math.sin(1.5 * self.time)))))
        star.move_to(50 * RIGHT)
        obs_points = Group(
            TrueDot(point, radius=0.1).set_color(GREEN).make_3d()
            for point in [orbit.get_top(), orbit.get_bottom()]
        )
        obs_lines = VGroup(
            self.get_obs_line(obs_point, star)
            for obs_point in obs_points
        )
        obs_lines.set_stroke(WHITE, 2)
        for line, point in zip(obs_lines, obs_points):
            line.start_point = point
            line.star = star
            line.add_updater(lambda m: m.put_start_and_end_on(m.start_point.get_center(), m.star.get_center()))

        obs_labels = VGroup(Text(f"Observação {n}") for n in [1, 2])
        for label, point, vect in zip(obs_labels, obs_points, [UP, DOWN]):
            label.next_to(point, vect, MED_SMALL_BUFF)

        self.add(star)

        self.play(
            ShowCreation(obs_lines[0], suspend_mobject_updating=True),
            FadeIn(obs_points[0]),
        )
        self.wait()
        self.play(Rotate(orbit, PI), run_time=2)
        self.play(
            ShowCreation(obs_lines[1], suspend_mobject_updating=True),
            FadeIn(obs_points[1]),
        )
        self.wait()

        # Show the angle vary during the orbit
        self.play(
            star.animate.move_to(15 * RIGHT),
            run_time=2
        )
        self.wait()

        


        curr_center = star.get_center()
        curr_angle = obs_lines[1].get_angle() - obs_lines[0].get_angle()
        # Label the distance and angle
        line_to_star = Line(sun.get_center(), star.get_center())
        line_to_star.add_updater(lambda m: m.put_start_and_end_on(sun.get_center(), star.get_center()))
        line_to_star.set_stroke(RED, 3)
        dist_label = Tex("D", font_size=60)
        dist_label.next_to(line_to_star, UP, buff=2 * SMALL_BUFF)
        dist_label.match_color(line_to_star)

        arc = Arc(PI, -curr_angle / 2, arc_center=star.get_center(), radius=3)
        arc.add_updater(lambda m: m.become(Arc(PI, -(obs_lines[1].get_angle() - obs_lines[0].get_angle()) / 2, arc_center=star.get_center(), radius=3)))
        arc_label = Tex(R"\theta / 2", font_size=60)
        arc_label.next_to(arc, LEFT, buff=SMALL_BUFF)

        self.play(
            ShowCreation(line_to_star),
            obs_lines.animate.set_stroke(width=1),
            FadeIn(dist_label, RIGHT),
        )
        self.wait()
        self.play(
            ShowCreation(arc),
        )
        self.wait()
        
        
        # # Pull it far away, then back
        curr_center = star.get_center()
        curr_angle = obs_lines[1].get_angle() - obs_lines[0].get_angle()
        orbit_radius / math.tan(curr_angle / 2)
        self.play(frame.animate.reorient(0,0,0,line_to_star.get_corner(RIGHT)))
        obs_lines.resume_updating()
        frame.add_updater(lambda m: m.reorient(0,0,0,line_to_star.get_corner(RIGHT)))
        self.play(
            UpdateFromAlphaFunc(star, lambda m, a: m.move_to(
                RIGHT * orbit_radius / math.tan(interpolate(curr_angle, 1e-5, there_and_back_with_pause(a)) / 2)
            )),
            
             run_time=15,
             rate_func=linear
        )
        
        
        # self.play(
        #     Transform(obs_lines[0].copy().clear_updaters(), obs_lines[1].copy(), remover=True),
        #     run_time=2
        # )
        # self.wait()

        #

    def get_obs_line(self, obj1, obj2, dash_length=0.1, stroke_color=WHITE, stroke_width=2):
        # line = DashedLine(obj1.get_center(), obj2.get_center())
        line = Line(obj1.get_center(), obj2.get_center())
        line.set_stroke(stroke_color, stroke_width)
        line.f_always.put_start_and_end_on(obj1.get_center, obj2.get_center)
        return line


class Raios23(Scene):
    def construct(self):
        self.frame.reorient(90,70,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3)
        esfera_celeste=EsferaCeleste()
        esfera_celeste.add_updater(lambda m: self.add(esfera_celeste))
        superficie = SuperficieObservador()
        self.add(superficie)
        self.add(esfera_celeste)
        sun = Sun(radius=RENDER_CELESTIAL_SPHERE_RADIUS*0.08).move_to(P(48,180).get_center())
        self.frame.add_ambient_rotation(4*DEG)
        self.play(FadeIn(sun))
        o = Line3D(ORIGIN,ORIGIN+[0,0,CELESTIAL_SPHERE_RADIUS*0.2],color=PURPLE,width=0.003)
        self.play(FadeIn(o))
        
        
        n_rays2 = 6
        rays = [Line3D(LEFT, RIGHT,width=0.03,color=YELLOW) for i in range(0,n_rays2)]
        group = [rays for i in range(0,n_rays2)]
        raios=Group()
        factor = 0.06
        for n, ray in enumerate(group):
            for i,r in enumerate(ray):
                r=Line3D(LEFT, RIGHT,width=0.01,color=YELLOW).set_opacity(0.7)
                r.put_start_and_end_on(
                    (n-(n_rays2-1)/2)*RIGHT*factor+(i-(n_rays2-1)/2)*DOWN*factor,
                    sun.get_center(),
                )
                raios.add(r)
        

        
        
        self.play(FadeIn(raios))
        self.wait(3)
        linha_astro = MarcadorAltura(P(48,180),arco=False,espessura_linha=5,cor_linha=RED)
        self.play(FadeIn(linha_astro))
        self.wait(5)
        raios.add_updater(lambda m: m.become(Group(*[
            Line3D(
                (n - (n_rays2 - 1) / 2) * RIGHT * factor + (i - (n_rays2 - 1) / 2) * DOWN * factor,
                sun.get_center(),
                width=0.001,
                color=YELLOW
            ).set_opacity(0.8)
            for n in range(n_rays2) for i in range(n_rays2)
        ])))
        self.play(sun.animate.shift(sun.get_center()*500 ),run_time=5)
        raios.clear_updaters()
        
        self.wait()
        sun2 = Sun(radius=RENDER_CELESTIAL_SPHERE_RADIUS*0.08).move_to(P(48,180).get_center())
        self.play(FadeIn(sun2))
        self.frame.add_ambient_rotation(-4*DEG)
        self.play(self.frame.animate.reorient(0,90-48,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3))
        self.wait(2)
        self.play(self.frame.animate.reorient(150,60,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3))
        self.wait(2)
        o2 = Line3D(ORIGIN,ORIGIN+[0,0,CELESTIAL_SPHERE_RADIUS*0.2],color=BLUE,width=0.003).shift(LEFT*0.12)
        self.play(FadeIn(o2))
        self.play(raios.animate.shift(RIGHT*0.12),o.animate.shift(RIGHT*0.12),o2.animate.shift(RIGHT*0.12))
        
        esfera_celeste2=EsferaCeleste().shift(RIGHT*0.12)
        self.play(FadeIn(esfera_celeste2))
        esfera_celeste2.add_updater(lambda m: self.add(esfera_celeste2))
        superficie2 = SuperficieObservador().shift(RIGHT*0.12)
        self.play(FadeIn(superficie2))
        sun2 = Sun(radius=RENDER_CELESTIAL_SPHERE_RADIUS*0.08).move_to(P(48,180).get_center()).shift(RIGHT*0.12)
        self.play(FadeIn(sun2))
        linha_astro2 = MarcadorAltura(P(48,180),arco=False,espessura_linha=5,cor_linha=RED).shift(RIGHT*0.12)
        self.play(FadeIn(linha_astro2))
        
        self.wait(2)
        stars_data = extract_star_data()

        stars = Group()

        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])/np.linalg.norm([x, y, z])*RENDER_CELESTIAL_SPHERE_RADIUS*500],
                color=color,
                radius=0.5+size/1200000,
                opacity=1
            )
            stars.add(star)
        stars.set_z_index(-3)
        esferaceleste = EsferaCeleste(opacidade=0.1,raio=RENDER_CELESTIAL_SPHERE_RADIUS*500)      
        self.play(FadeIn(esferaceleste),FadeIn(stars))
        esferaceleste.add_updater(lambda m: self.add(esferaceleste))
        
        sun.scale(200)
        earth=Earth(radius=RENDER_CELESTIAL_SPHERE_RADIUS*30,clouds=False).rotate(110*DEG,UP).rotate(170*DEG,OUT).rotate(30*DEG,RIGHT).next_to(superficie.get_center(),IN*RENDER_CELESTIAL_SPHERE_RADIUS*30,buff=0)
        self.play(FadeIn(earth))
        self.wait(2)
        self.play(self.frame.animate.reorient(150,60,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1200),run_time=10,rate_func=smooth)
        self.frame.add_ambient_rotation(4*DEG)
        self.wait(4)