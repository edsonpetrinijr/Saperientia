from manimlib import *
from manimlib.Saperientia.astro_structures import * 
from manimlib.Saperientia.stellarium import * 
from manimlib.Saperientia.Variables import * 
import numpy as np
from scipy.spatial.transform import Rotation as R

class TimeConversion(Scene):
    def construct(self):
        # First conversion: 1h -> 60min
        hour_text = Tex("1h").scale(2).to_edge(UL).shift(RIGHT)
        minutes_text_60 = Tex("60min").scale(2).next_to(hour_text, DOWN, buff=LARGE_BUFF)

        self.play(Write(hour_text))
        self.wait()

        # Transform a copy of "1h" to "60min"
        self.play(TransformFromCopy(hour_text, minutes_text_60))
        self.wait()

        # Second conversion: 1min -> 60s
        minutes_text_1 = Tex("1min").scale(2).to_edge(UP).shift(RIGHT)
        seconds_text_60 = Tex("60s").scale(2).next_to(minutes_text_1, DOWN, buff=LARGE_BUFF)


        self.play(Write(minutes_text_1))
        self.wait()

        # Transform a copy of "1min" to "60s"
        self.play(TransformFromCopy(minutes_text_1, seconds_text_60))
        self.wait(10)

class DemonstrateAnArcSecond(InteractiveScene):
    def construct(self):
        # Add circle
        frame = self.frame
        frame.set_height(90)
        frame.set_field_of_view(10 * DEG)
        radius = 35

        circle = Circle(radius=radius)
        circle.set_stroke(WHITE, 2)

        frame = self.frame
        deg_height = radius * DEG
        zero_point = circle.get_right()

        # Add degree tick marks
        deg_ticks = VGroup(
            self.get_tick_mark(circle, n * DEG, 1, max_aparent_length=1)
            for n in range(360)
        )
        deg10_ticks = VGroup(
            self.get_tick_mark(circle, n * DEG, 2.5, max_aparent_length=1)
            for n in range(0, 360, 10)
        )
        deg30_labels = VGroup(
            self.get_tick_label(tick, number)
            for tick, number in zip(deg10_ticks[0::3], range(0, 360, 30))
        )
        deg10_labels = VGroup(
            self.get_tick_label(tick, number, visibility_range=(60, 30))
            for tick, number in zip(
                [*deg10_ticks[-2:], *deg10_ticks[1:3]],
                [340, 350, 10, 20]
            )
        )
        deg_labels = VGroup(
            self.get_tick_label(tick, number, visibility_range=(12, 6))
            for tick, number in zip(deg_ticks[1:10], range(1, 10))
        )

        self.add(circle, deg_ticks, deg10_ticks)
        self.add(deg30_labels, deg10_labels, deg_labels)

        # Add arc minute labels
        arc_minute_ticks = VGroup(
            self.get_tick_mark(circle, n * DEG / 60, 2.5e-2, max_aparent_length=0.75)
            for n in range(0, 60)
        )
        arc_minute_labels = VGroup(
            self.get_tick_label(
                tick,
                number,
                unit="minutos de arco",
                frame_prop=0.025,
                visibility_range=(0.35 * deg_height, 0.15 * deg_height),
            )
            for tick, number in zip(arc_minute_ticks[1:20], range(1, 20))
        )

        self.add(arc_minute_ticks)
        self.add(arc_minute_labels)

        # Add arc second labels
        arc_second_ticks = VGroup(
            self.get_tick_mark(circle, n * DEG / 60 / 60, 3.5e-4, max_aparent_length=0.5)
            for n in range(0, 60)
        )
        arc_second_ticks.rotate(0.25 * DEG, about_point=zero_point)

        am_height = deg_height / 60
        arc_second_labels = VGroup(
            self.get_tick_label(
                tick,
                number,
                unit="segundos de arco",
                frame_prop=0.025,
                visibility_range=(0.35 * am_height, 0.15 * am_height),
                font_size=0.1
            )
            for tick, number in zip(arc_second_ticks[1:20], range(1, 20))
        )

        self.add(arc_second_ticks)
        self.add(arc_second_labels)

        # Zoom in to one degree


        self.play(
            frame.animate.set_height(1.5 * deg_height).move_to(circle.pfp(0.5 * DEG / TAU)),
            run_time=8
        )
        self.wait(2)



        # Zoom in to arc seconds
        self.play(
            frame.animate.set_height(1.5 * deg_height / 60).move_to(circle.pfp(0.5 * DEG / TAU / 60)),
            run_time=4
        )
        self.wait()

        op_tracker = ValueTracker(0)
        arc_second_labels.add_updater(lambda m: m.set_opacity(op_tracker.get_value()))
        self.play(
            VGroup(circle, arc_minute_ticks, arc_second_ticks).animate.scale(20, about_point=zero_point),
            op_tracker.animate.set_value(1).set_anim_args(time_span=(1, 3)),
            run_time=4
        )
        self.wait()

    def get_tick_mark(
        self,
        circle,
        angle,
        tick_length,
        max_aparent_length=0.5,
        stroke_color=WHITE,
        stroke_width=2
    ):
        tick = Line(ORIGIN, tick_length * RIGHT)
        tick.move_to(circle.get_right())
        tick.rotate(angle, about_point=circle.get_center())
        tick.set_stroke(WHITE, stroke_width)

        frame_prop = max_aparent_length / FRAME_WIDTH
        tick.add_updater(lambda m: m.set_length(
            min(frame_prop * self.frame.get_width(), tick_length)
        ))

        return tick

    def get_tick_label(
        self,
        tick,
        number,
        unit=R"^\circ",
        buff_ratio=0.5,
        font_size=360,
        visibility_range=(1000, 900),
        frame_prop=0.05,
    ):
        if not unit.startswith("^"):
            # label[-1].shift(0.5 * label[0].get_width() * RIGHT)
            if number == 1:
                if unit in "minutos de arco":
                    unit = "minuto de arco"
                elif unit in "segundos de arco":
                    unit = "segundo de arco"
        label = Integer(number, unit=unit, font_size=font_size)

        direction = normalize(tick.get_vector())
        direction[abs(direction) < 0.2] = 0

        label_height = label.get_height()

        def update_label(label):
            frame_height = self.frame.get_height()

            # Opacity
            opacity = float(clip(inverse_interpolate(*visibility_range, frame_height), 0, 1))
            label.set_opacity(opacity)

            # Max height
            height = min(frame_prop * frame_height, label_height)
            label.set_height(height)

            # Location
            buff = buff_ratio * label[0].get_width()
            label.next_to(tick.get_end(), direction, buff=buff)
            return label

        label.add_updater(update_label)

        return label

class ArcSecondConversion(Scene):
    def construct(self):
        # Starting with 20 arc-seconds
        tex2 = Tex(r"10'").scale(2).to_edge(UP)
        tex1 = Tex(r"5^{\circ}").scale(2).next_to(tex2,LEFT,buff=LARGE_BUFF)
        tex3 = Tex(r"15''").scale(2).next_to(tex2,RIGHT,buff=LARGE_BUFF)
        
        self.play(Write(tex2),Write(tex1),Write(tex3))
        self.wait(2)
        self.play(tex3.animate.move_to([0,0,0]))
        
        sum_line = Line(
            tex3.get_corner(UR) + RIGHT*0.2, 
            tex3.get_corner(DR) + RIGHT*0.2+DOWN*0.2
        )
        sum_line2 = Line(
            tex3.get_corner(DR) + RIGHT*3+DOWN*0.2, 
            tex3.get_corner(DR) + RIGHT*0.2+DOWN*0.2
        )
        self.play(ShowCreation(sum_line),ShowCreation(sum_line2))
        sixty = Tex(r"60").scale(2).next_to(sum_line2,UP)
        self.play(Write(sixty))
        self.wait(2)
        quarter = Tex(r"0,25'").scale(2).next_to(sum_line2,DOWN)
        self.play(Write(quarter))
        self.wait(2)
        self.play(FadeOut(sum_line),FadeOut(sum_line2),FadeOut(sixty),FadeOut(tex3))
        self.play(quarter.animate.move_to(tex2.get_corner(DOWN)+DOWN*1))
        plus = Tex(r"+").scale(2).next_to(quarter,LEFT)
        self.play(Write(plus))
        sum_line3 = Line(
            quarter.get_corner(DR) + RIGHT*0.2, 
            quarter.get_corner(DL) + LEFT*0.2
        )
        self.play(ShowCreation(sum_line3))
        self.wait(2)
        newminute = Tex(r"10,25'").scale(2).next_to(sum_line3,DOWN)
        self.play(Write(newminute))
        self.wait(2)
        self.play(FadeOut(tex2), FadeOut(plus),FadeOut(sum_line3),FadeOut(quarter))
        
        sum_line4 = Line(
            newminute.get_corner(UR) + RIGHT*0.2, 
            newminute.get_corner(DR) + RIGHT*0.2+DOWN*0.2
        )
        sum_line5 = Line(
            newminute.get_corner(DR) + RIGHT*3+DOWN*0.2, 
            newminute.get_corner(DR) + RIGHT*0.2+DOWN*0.2
        )
        self.play(ShowCreation(sum_line4),ShowCreation(sum_line5))
        sixty = Tex(r"60").scale(2).next_to(sum_line5,UP)
        self.play(Write(sixty))
        self.wait(2)
        newdegree = Tex(r"0,17083^{\circ}").scale(2).next_to(sum_line5,DOWN)
        self.play(Write(newdegree))
        self.wait(2)
        self.play(FadeOut(sum_line4),FadeOut(sum_line5),FadeOut(sixty),FadeOut(newminute))
        self.play(newdegree.animate.move_to(tex1.get_corner(DOWN)+DOWN*1))
        plus = Tex(r"+").scale(2).next_to(newdegree,LEFT)
        self.play(Write(plus))
        sum_line6 = Line(
            newdegree.get_corner(DR) + RIGHT*0.2, 
            newdegree.get_corner(DL) + LEFT*0.2
        )
        self.play(ShowCreation(sum_line6))
        self.wait(2)
        newatall = Tex(r"5,17083^{\circ}").scale(2).next_to(sum_line6,DOWN)
        self.play(Write(newatall))
        self.play(FadeOut(newdegree),FadeOut(plus),FadeOut(sum_line6),FadeOut(tex1))
        self.wait(5)
        
class Conversion2(Scene):
    def construct(self):
        # Começando com 30 segundos de arco
        arcseconds = Tex(r"30''").scale(2).shift(LEFT*6)
        self.play(Write(arcseconds))
        self.wait(1)
        
        # Primeira seta curvada (segundos para minutos)
        arrow1 = CurvedArrow(
            arcseconds.get_right() + RIGHT*0.5,
            arcseconds.get_right() + RIGHT*2.5,
            angle=-PI/3
        )
        
        # Texto "÷60" embaixo da primeira seta
        div60_1 = Text("÷60").scale(1)
        div60_1.next_to(arrow1, DOWN, buff=0.3)
        
        # Resultado em minutos
        minutes = Tex(r"0,5'").scale(2)
        minutes.move_to(arcseconds.get_center() + RIGHT*4.5)
        
        self.play(ShowCreation(arrow1))
        self.play(Write(div60_1))
        self.wait(1)
        self.play(Write(minutes))
        self.wait(1)
        
        # Segunda seta curvada (minutos para graus)
        arrow2 = CurvedArrow(
            minutes.get_right() + RIGHT*0.5,
            minutes.get_right() + RIGHT*2.5,
            angle=-PI/3
        )
        
        # Texto "÷60" embaixo da segunda seta
        div60_2 = Text("÷60").scale(1)
        div60_2.next_to(arrow2, DOWN, buff=0.3)
        
        # Resultado em graus
        degrees = Tex(r"0,00833^{\circ}").scale(2)
        degrees.move_to(minutes.get_center() + RIGHT*5.5)
        
        self.play(ShowCreation(arrow2))
        self.play(Write(div60_2))
        self.wait(1)
        self.play(Write(degrees))
        self.wait(3)

class ConversionTime(Scene):
    def construct(self):
        # Começando com 30 segundos de arco
        arcseconds = Tex(r"30s").scale(2).shift(LEFT*6)
        self.play(Write(arcseconds))
        self.wait(1)
        
        # Primeira seta curvada (segundos para minutos)
        arrow1 = CurvedArrow(
            arcseconds.get_right() + RIGHT*0.5,
            arcseconds.get_right() + RIGHT*2.5,
            angle=-PI/3
        )
        
        # Texto "÷60" embaixo da primeira seta
        div60_1 = Text("÷60").scale(1)
        div60_1.next_to(arrow1, DOWN, buff=0.3)
        
        # Resultado em minutos
        minutes = Tex(r"0,5min").scale(2)
        minutes.move_to(arcseconds.get_center() + RIGHT*5.5)
        
        self.play(ShowCreation(arrow1))
        self.play(Write(div60_1))
        self.wait(1)
        self.play(Write(minutes))
        self.wait(1)
      
class ConversionTime2(Scene):
    def construct(self):
        # Começando com 30 segundos de arco
        arcseconds = Tex(r"30min").scale(2).shift(LEFT*4)
        self.play(Write(arcseconds))
        self.wait(1)
        
        # Primeira seta curvada (segundos para minutos)
        arrow1 = CurvedArrow(
            arcseconds.get_right() + RIGHT*0.5,
            arcseconds.get_right() + RIGHT*2.5,
            angle=-PI/3
        )
        
        # Texto "÷60" embaixo da primeira seta
        div60_1 = Text("÷60").scale(1)
        div60_1.next_to(arrow1, DOWN, buff=0.3)
        
        # Resultado em minutos
        minutes = Tex(r"0,5h").scale(2)
        minutes.move_to(arcseconds.get_center() + RIGHT*5.5)
        
        self.play(ShowCreation(arrow1))
        self.play(Write(div60_1))
        self.wait(1)
        self.play(Write(minutes))
        self.wait(1)
      
class DecimalConversion(Scene):
    def construct(self):
        # Começando com 5,17083 graus
        degrees = Tex(r"5",r"\phantom{000000}",r"^{\circ}").scale(2).shift(LEFT*4+UP)
        degrees2= Tex(r",17083").scale(2).move_to(degrees.get_center()+DOWN*0.1)
        self.play(Write(degrees), Write(degrees2))
        self.wait(1)
        
        # Separar parte inteira e decimal
        whole_part = Tex(r"5",r"^{\circ}").scale(2)
        whole_part.move_to(degrees.get_center() + UP*2)
        
        decimal_part = Tex(r"0",r",17083").scale(2)
        decimal_part.move_to(degrees.get_center() + DOWN*1)
        
        self.play(TransformMatchingTex(degrees,whole_part),TransformMatchingTex(degrees2,decimal_part))
        self.wait(1)
        
        # Primeira seta curvada (decimal × 60)
        arrow1 = CurvedArrow(
            decimal_part.get_right() + RIGHT*0.5,
            decimal_part.get_right() + RIGHT*2.5,
            angle=-PI/3
        )
        
        # Texto "×60" embaixo da primeira seta
        mult60_1 = Tex(r"\times 60").scale(1)
        mult60_1.next_to(arrow1, DOWN, buff=0.3)
        
        # Resultado em minutos (0,17083 × 60 = 10,2498)
        minutes_result = Tex(r"10",r"\phantom{1111}",r"'").scale(2)
        minutes_result.move_to(decimal_part.get_center() + RIGHT*5.8)
        minutes_result2 = Tex(r",25").scale(2).move_to(minutes_result.get_center()+DOWN*0.1)
        
        self.play(ShowCreation(arrow1))
        self.play(Write(mult60_1))
        self.wait(1)
        self.play(Write(minutes_result),Write(minutes_result2))
        self.wait(1)
        
        # Separar parte inteira e decimal dos minutos
        minutes_whole = Tex(r"10",r"'").scale(2)
        minutes_whole.move_to(minutes_result.get_center() + UP*1.5)
        
        minutes_decimal = Tex(r"0",r",25").scale(2)
        minutes_decimal.move_to(minutes_result.get_center() + DOWN*1)
        
        self.play(TransformMatchingTex(minutes_result,minutes_whole),TransformMatchingTex(minutes_result2,minutes_decimal))
        self.wait(1)
        
        # Segunda seta curvada (decimal × 60)
        arrow2 = CurvedArrow(
            minutes_decimal.get_right() + RIGHT*0.5,
            minutes_decimal.get_right() + RIGHT*2.5,
            angle=-PI/3
        )
        
        # Texto "×60" embaixo da segunda seta
        mult60_2 = Tex(r"\times 60").scale(1)
        mult60_2.next_to(arrow2, DOWN, buff=0.3)
        
        # Resultado em segundos (0,2498 × 60 = 14,988)
        seconds_result = Tex(r"15''").scale(2)
        seconds_result.move_to(minutes_decimal.get_center() + RIGHT*4.2)
        
        self.play(ShowCreation(arrow2))
        self.play(Write(mult60_2))
        self.wait(1)
        self.play(Write(seconds_result))
        self.wait(1)
        
        # Resultado final
        final_result = VGroup(
            Tex(r"5^{\circ}").scale(2.5),
            Tex(r"10'").scale(2.5),
            Tex(r"15''").scale(2.5)
        ).arrange(RIGHT, buff=0.5).move_to(DOWN*3)
        
        self.play(Write(final_result))
        self.wait(3)
                
class Hmins(Scene):
    def construct(self):
        final_result = VGroup(
            Tex(r"h ").scale(2.5),
            Tex(r"min").scale(2.5),
            Tex(r"s").scale(2.5)
        ).arrange(RIGHT, buff=1).move_to(DOWN*3)
        
        self.play(Write(final_result))
        self.wait(3)
        
class UnitsAngle(Scene):
    def construct(self):
        final_result = VGroup(
            Tex(r"^{\circ}").scale(2.5),
            Tex(r"'").scale(2.5),
            Tex(r"''").scale(2.5)
        ).arrange(RIGHT, buff=1).move_to(DOWN*3)
        
        self.play(Write(final_result))
        self.wait(3)

class VerticalConversion(Scene):
    def construct(self):
        # 1 hora = 15 graus (no topo)
        hour_text = Tex(r"1h").scale(2).move_to(UP*3 + LEFT*2)
        self.play(Write(hour_text))
        self.wait(1)
        
        # Seta para a direita
        arrow1 = Arrow(
            hour_text.get_right() + RIGHT*0.2,
            hour_text.get_right() + RIGHT*2,
            buff=0
        )
        self.play(ShowCreation(arrow1))
        
        # 15 graus
        degrees_text = Tex(r"15^{\circ}").scale(2).move_to(UP*3 + RIGHT*2)
        self.play(Write(degrees_text))
        self.wait(1)
        
        # 1 minuto = 15 minutos de arco (no meio)
        minute_text = Tex(r"1min").scale(2).move_to(UP*1 + LEFT*2)
        self.play(Write(minute_text))
        self.wait(1)
        
        # Seta para a direita
        arrow2 = Arrow(
            minute_text.get_right() + RIGHT*0.2,
            minute_text.get_right() + RIGHT*2,
            buff=0
        )
        self.play(ShowCreation(arrow2))
        
        # 15 minutos de arco
        arcminutes_text = Tex(r"15'").scale(2).move_to(UP*1 + RIGHT*2)
        self.play(Write(arcminutes_text))
        self.wait(1)
        
        # 1 segundo = 15 segundos de arco (embaixo)
        second_text = Tex(r"1s").scale(2).move_to(DOWN*1 + LEFT*2)
        self.play(Write(second_text))
        self.wait(1)
        
        # Seta para a direita
        arrow3 = Arrow(
            second_text.get_right() + RIGHT*0.2,
            second_text.get_right() + RIGHT*2,
            buff=0
        )
        self.play(ShowCreation(arrow3))
        
        # 15 segundos de arco
        arcseconds_text = Tex(r"15''").scale(2).move_to(DOWN*1 + RIGHT*2)
        self.play(Write(arcseconds_text))
        self.wait(3) 
    
class Episodio1(Scene):
    def construct(self):
        self.frame.reorient(0, 0, 0, ORIGIN+[0,0,0.001], 0.0005)
        # self.camera.frame.set_width(10000)
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
    
        sun = Sun()
        sun.add_updater(lambda m:self.add(sun))
        earth = Earth(clouds=False).move_to(RENDER_SUN_EARTH_DISTANCE*RIGHT)
        light = self.camera.light_source
        moon = Moon()  # Lua ainda menor
        mars = Mars()  # Lua ainda menor

        self.add(sun,mars)
        
        

        # Trackers para os ângulos da Terra e da Lua
        angle_earth = ValueTracker(0)
        angle_moon = ValueTracker(0)

        # Atualiza a posição da Terra orbitando o Sol
        # earth_orbit = always_redraw(lambda: 
        #     earth.move_to(sun.get_center() + RENDER_SUN_EARTH_DISTANCE * np.array([
        #         np.cos(angle_earth.get_value()),
        #         np.sin(angle_earth.get_value()),
        #         0
        #     ]))
        # )
        mars_orbit = always_redraw(lambda: 
            mars.move_to(sun.get_center() + RENDER_SUN_MARS_DISTANCE * np.array([
                np.cos(angle_earth.get_value()),
                np.sin(angle_earth.get_value()),
                0
            ]))
        )
        
        # Atualiza a posição da Lua orbitando a Terra
        mars_orbit = always_redraw(lambda: 
            mars.move_to(sun.get_center() + RENDER_SUN_MARS_DISTANCE * np.array([
                np.cos(angle_earth.get_value()),
                np.sin(angle_earth.get_value()) * np.cos(np.radians(23.5)),
                np.sin(angle_earth.get_value()) * np.sin(np.radians(23.5))
            ]))
        )
        # self.frame.reorient(330,70,0,[53.38888931,  0.        ,  0.        ],15)
        self.frame.reorient(340,105,0,earth.get_center(),1)
        self.wait(5)
        self.play(
            angle_earth.animate.increment_value( PI),
            angle_moon.animate.increment_value(8 * PI),
            # Rotate(earth,10*TAU,axis=Z_AXIS),
            run_time=20,
            rate_func=linear
        )
 
class CenaEsfera(Scene):
    def construct(self):
        self.frame.reorient(0, 80, 0, ORIGIN+[0,0,0.001], 0.05)
        # self.camera.frame.set_width(10000)
        stars_data = extract_star_data()

        stars = Group()

        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])/np.linalg.norm([x, y, z])*RENDER_CELESTIAL_SPHERE_RADIUS],
                color=color,
                radius=0.0007+size/800000000,
                opacity=1
            )
            stars.add(star)

        self.add(stars)
        esferaceleste = EsferaCeleste(opacidade=0.1)        
        sun = Sun(radius=0.07*ELEMENTS_SCALE, big_glow_ratio=0.1*RENDER_SUN_RADIUS).move_to(esferaceleste.get_corner(RIGHT))
        sun.add_updater(lambda m:self.add(sun))
        self.add(esferaceleste)
        esferaceleste.add_updater(lambda m:self.add(esferaceleste))
        earth = Earth(clouds=False).scale(0.2)
        light = self.camera.light_source
        moon = Moon().scale(0.1).move_to(esferaceleste.get_corner(LEFT)) # Lua ainda menor
        mars = Mars().scale(0.1).move_to(esferaceleste.get_corner(UP))  # Lua ainda menor
        self.add(sun,earth,moon,mars)
        self.play(self.frame.animate.reorient(180, 80, 0, ORIGIN, 0.05),run_time=5)
        self.play(self.frame.animate.reorient(90, 80, 0, ORIGIN, 0.5),run_time=5)

class Escala(Scene):
    def construct(self):
        self.frame.reorient(0, 0, 0, ORIGIN+[0,0,0.001], 0.0005)
        # self.camera.frame.set_width(10000)
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
    
        sun = Sun()
        sun.add_updater(lambda m:self.add(sun))
        earth = Earth(clouds=False).move_to(RENDER_SUN_EARTH_DISTANCE*RIGHT)
        light = self.camera.light_source
        moon = Moon()  # Lua ainda menor
        mars = Mars()  # Lua ainda menor

        self.add(sun,mars,earth,moon)
        
        

        # Trackers para os ângulos da Terra e da Lua
        angle_earth = ValueTracker(0)
        angle_moon = ValueTracker(0)


        earth_orbit = always_redraw(lambda: 
            earth.move_to(sun.get_center() + RENDER_SUN_EARTH_DISTANCE * np.array([
                np.cos(angle_earth.get_value()),
                np.sin(angle_earth.get_value()),
                0
            ]))
        )
        
        # Atualiza a posição da Lua orbitando a Terra
        moon_orbit = always_redraw(lambda:
            moon.move_to(earth.get_center() + RENDER_EARTH_MOON_DISTANCE * np.array([
                np.cos(angle_moon.get_value()),
                np.sin(angle_moon.get_value()),
                0
            ]))
        )
        # Atualiza a posição da Lua orbitando a Terra
        mars_orbit = always_redraw(lambda: 
            mars.move_to(sun.get_center() + RENDER_SUN_MARS_DISTANCE * np.array([
                np.cos(angle_earth.get_value()),
                np.sin(angle_earth.get_value()) * np.cos(np.radians(23.5)),
                np.sin(angle_earth.get_value()) * np.sin(np.radians(23.5))
            ]))
        )
        # self.frame.reorient(330,70,0,[53.38888931,  0.        ,  0.        ],15)
        self.frame.reorient(340,105,0,earth.get_center(),1)
        self.play(self.frame.animate.reorient(270,105,0,earth.get_center(),1),run_time=5)
        self.play(self.frame.animate.reorient(270, 80, 0, moon.get_center(), 0.5),run_time=3)
        self.play(self.frame.animate.reorient(90, 80, 0, moon.get_center(), 0.5),run_time=2)
        self.play(self.frame.animate.reorient(90, 80, 0, sun.get_center(), 20),run_time=2)
        self.wait(1)
        self.play(self.frame.animate.reorient(270, 80, 0, mars.get_center(), 1),run_time=2)
        self.wait(2)
        self.play(self.frame.animate.reorient(0, 90, 0, mars.get_center(), 1),run_time=2)
        self.play(self.frame.animate.reorient(0, 91, 0, [627751.2200000001, 6026669.100000001, -127126.83], 1000000),run_time=30,rate_func=smooth)

class CenaEsfera2(Scene):
    def construct(self):
        self.frame.reorient(0, 80, 0, ORIGIN+[0,0,0.001], 0.05)
        # self.camera.frame.set_width(10000)
        stars_data = extract_star_data()

        stars = Group()

        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])/np.linalg.norm([x, y, z])*RENDER_CELESTIAL_SPHERE_RADIUS],
                color=color,
                radius=0.0007+size/800000000,
                opacity=1
            )
            stars.add(star)

        self.add(stars)
        esferaceleste = EsferaCeleste(opacidade=0.1)        
        sun = Sun(radius=0.07*ELEMENTS_SCALE, big_glow_ratio=0.1*RENDER_SUN_RADIUS).move_to(esferaceleste.get_corner(RIGHT))
        sun.add_updater(lambda m:self.add(sun))
        self.add(esferaceleste)
        esferaceleste.add_updater(lambda m:self.add(esferaceleste))
        earth = Earth(clouds=False).scale(0.2)
        light = self.camera.light_source
        moon = Moon().scale(0.1).move_to(esferaceleste.get_corner(LEFT)) # Lua ainda menor
        mars = Mars().scale(0.1).move_to(esferaceleste.get_corner(UP))  # Lua ainda menor
        self.add(sun,earth,moon,mars)
        self.play(self.frame.animate.reorient(180, 80, 0, ORIGIN, 0.05),run_time=5)
        self.play(self.frame.animate.reorient(120, 50, 0, ORIGIN, 0.03),run_time=3)
        self.play(self.frame.animate.reorient(-60, 110, 0, ORIGIN, 0.03),run_time=5)
        self.play(self.frame.animate.reorient(180, 80, 0, ORIGIN, 0.05),run_time=5)
        self.play(earth.animate.scale(0.01),run_time=2)
        self.play(self.frame.animate.reorient(90, 80, 0, ORIGIN, 0.5),run_time=5)
        self.play(self.frame.animate.reorient(180, 80, 0, ORIGIN, 0.5),run_time=5)
        self.play(self.frame.animate.reorient(180, 80, 0, ORIGIN, 0.05),run_time=5)
        self.wait(4)