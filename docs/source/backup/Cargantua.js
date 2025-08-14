(async () => {
    try {
        const { default: canvasManager } = await import('./canvasmanager.js');
        const { ShaderBuilder } = await import('./shaderbuilder.js');
        const { default: FPSCounter } = await import('./fps.js');

        new FPSCounter();

        const manager = canvasManager().attach('shader');
        const gl = manager.context('webgl2');
        manager.resize('full', 'full');
        manager.listen('resize', 250).listen('dpr', 250);

        const frag = `#version 300 es
        precision highp float;

        uniform vec3 u_resolution;
        uniform float u_time;
        uniform float u_mix;
        out vec4 fragColor;

        void main() {
            vec2 fragCoord = gl_FragCoord.xy;
            vec2 r = u_resolution.xy;
            vec2 p = (fragCoord + fragCoord - r) / r.y;
            vec2 z = vec2(0.5);
            vec2 i = vec2(0.1);
            vec2 f = p * (z += 5. - 6. * exp(.4 - dot(p, p)));
            vec4 O = vec4(0.0);
            for (i.y = 1.0; i.y <= 8.0; i.y += 1.0) {
                O += (tanh(f) + 1.0).xyyx * abs(f.x - f.y);
                f += tanh(f.yx * i.y + i + u_time) / i.y + 0.7;
            }
            O = tanh(5.0 * exp(z.x - 4.0 - p.y * vec4(-1.0, 1.0, 2.0, 0.0)) / O);

            float mixPhase = dot(p, p) + z.x + u_time + sin(p.x * 1.5 + p.y * 2.5 + u_time * 0.5);
            float channel = cos(mixPhase * 4.0);
            vec3 glow = vec3(
                0.6 + 0.4 * sin(channel + 1.0),
                0.6 + 0.4 * sin(channel + 0.0),
                0.6 + 0.4 * sin(channel + 2.0)
            );
            O.rgb *= glow * 1.;

            fragColor = O;
        }`;

        const vert = `#version 300 es
        precision highp float;

        layout(location = 0) in vec2 a_position;
        out vec2 vUv;
        void main() {
            gl_Position = vec4(a_position, 0.0, 1.0);
            vUv = a_position * 0.5 + 0.5;
        }`;

        const shader = new ShaderBuilder(gl, vert, frag);
        shader.build();

        const quadVerts = new Float32Array([
            -1, -1,
            1, -1,
            -1,  1,
            1,  1,
        ]);
        shader.createVAO('quad');
        shader.bindVAO('quad');
        shader.createVBO('quad', quadVerts);
        shader.setAttribute('a_position', 2, gl.FLOAT);

        const draw = (time) => {
            // Clear the canvas
            gl.clearColor(0, 0, 0, 1);
            gl.clear(gl.COLOR_BUFFER_BIT);

            shader.use();
            shader.setUniform3f('u_resolution', gl.canvas.width, gl.canvas.height, 1);
            shader.setUniform1f('u_time', time * 0.001);
            shader.setUniform1f('u_mix', 0.5);

            gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

            requestAnimationFrame(draw);
        };

        requestAnimationFrame(draw);

    } catch (error) {
        console.error('Cargantua.js error:', error);
        // Fallback: create a simple color animation
        const canvas = document.getElementById('shader');
        if (canvas) {
            const ctx = canvas.getContext('2d');
            canvas.width = window.innerWidth;
            canvas.height = window.innerHeight;

            const animate = (time) => {
                const hue = (time * 0.05) % 360;
                ctx.fillStyle = `hsl(${hue}, 50%, 10%)`;
                ctx.fillRect(0, 0, canvas.width, canvas.height);
                requestAnimationFrame(animate);
            };
            requestAnimationFrame(animate);
        }
    }
})();
