/**
 * ShaderBuilder: Enhanced utility for compiling, linking, and managing WebGL shaders with multi-pass support
 * Provides uniform/attribute setting, buffer binding, logging, validation, resource cleanup, and pipeline management.
 */
const ShaderBuilder = class {
    constructor(gl, vsource, fsource) {
        if (!gl || !(gl instanceof WebGL2RenderingContext)) {
            throw new Error("ShaderBuilder requires a WebGL2 context.");
        }
        if (!vsource.includes("#version 300 es") || !fsource.includes("#version 300 es")) {
            throw new Error("Shader sources must include '#version 300 es'.");
        }
        this.gl = gl;
        this.vsource = vsource;
        this.fsource = fsource;
        this.program = null;
        this.debug = false;
        this._uniformCache = {};
        this._attributeCache = {};
        this._ubos = new Map();
        this._vaos = new Map();
        this._vbos = new Map();
        this._textures = new Map();
        this._framebuffers = new Map();
        this._activeTextureUnit = 0;
    }

    setDebug(debug) {
        if (typeof debug !== 'boolean') throw new Error("Debug flag must be a boolean.");
        this.debug = debug;
    }

    setVertexSource(source) {
        if (typeof source !== "string" || source.trim() === "") throw new Error("Vertex shader source must be a non-empty string.");
        this.vsource = source;
        if (this.debug) console.debug('[ShaderBuilder] setVertexSource called.', source);
    }

    setFragmentSource(source) {
        if (typeof source !== "string" || source.trim() === "") throw new Error("Fragment shader source must be a non-empty string.");
        this.fsource = source;
        if (this.debug) console.debug('[ShaderBuilder] setFragmentSource called.', source);
    }

    build(debug = false) {
        const gl = this.gl;
        if (!this.vsource || !this.fsource) throw new Error("Shader sources must not be empty.");

        // Clean up existing program if rebuilding
        if (this.program) {
            gl.deleteProgram(this.program);
        }

        const vertexShader = this._compileShader(gl.VERTEX_SHADER, this.vsource);
        const fragmentShader = this._compileShader(gl.FRAGMENT_SHADER, this.fsource);
        this.program = gl.createProgram();
        gl.attachShader(this.program, vertexShader);
        gl.attachShader(this.program, fragmentShader);
        gl.linkProgram(this.program);
        gl.deleteShader(vertexShader);
        gl.deleteShader(fragmentShader);
        if (!gl.getProgramParameter(this.program, gl.LINK_STATUS)) {
            const info = gl.getProgramInfoLog(this.program);
            console.error("[ShaderBuilder] Failed to link shader program:", info);
            throw new Error("Failed to link shader program: " + info);
        }
        this._uniformCache = {};
        this._attributeCache = {};
        this._ubos.clear();

        if (this.debug) console.debug('[ShaderBuilder] Shader program built successfully.');
    }

    // Hot reloading method
    rebuild() {
        if (this.debug) console.debug('[ShaderBuilder] Rebuilding shader program...');
        this.build();
        return this;
    }

    // Update both shaders and rebuild
    updateShaders(vsource, fsource) {
        if (vsource) this.setVertexSource(vsource);
        if (fsource) this.setFragmentSource(fsource);
        return this.rebuild();
    }

    _compileShader(type, source) {
        const gl = this.gl;
        const shader = gl.createShader(type);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            const info = gl.getShaderInfoLog(shader);
            console.error(`[ShaderBuilder] Failed to compile shader:\n${source}\nError: ${info}`);
            throw new Error("Failed to compile shader: " + info);
        }
        return shader;
    }

    // Uniform helper methods
    setUniform1f(name, value) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform1f(location, value);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform1f: ${name} = ${value}`);
        }
    }

    setUniform2f(name, x, y) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform2f(location, x, y);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform2f: ${name} = (${x}, ${y})`);
        }
    }

    setUniform3f(name, x, y, z) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform3f(location, x, y, z);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform3f: ${name} = (${x}, ${y}, ${z})`);
        }
    }

    setUniform4f(name, x, y, z, w) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform4f(location, x, y, z, w);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform4f: ${name} = (${x}, ${y}, ${z}, ${w})`);
        }
    }

    setUniform1i(name, value) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform1i(location, value);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform1i: ${name} = ${value}`);
        }
    }

    setUniform2i(name, x, y) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform2i(location, x, y);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform2i: ${name} = (${x}, ${y})`);
        }
    }

    setUniform3i(name, x, y, z) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform3i(location, x, y, z);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform3i: ${name} = (${x}, ${y}, ${z})`);
        }
    }

    setUniform4i(name, x, y, z, w) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform4i(location, x, y, z, w);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform4i: ${name} = (${x}, ${y}, ${z}, ${w})`);
        }
    }

    setUniformMatrix2fv(name, matrix) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniformMatrix2fv(location, false, matrix);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniformMatrix2fv: ${name}`, matrix);
        }
    }

    setUniformMatrix3fv(name, matrix) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniformMatrix3fv(location, false, matrix);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniformMatrix3fv: ${name}`, matrix);
        }
    }

    setUniformMatrix4fv(name, matrix) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniformMatrix4fv(location, false, matrix);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniformMatrix4fv: ${name}`, matrix);
        }
    }

    setUniformArray1fv(name, array) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform1fv(location, array);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform1fv: ${name}`, array);
        }
    }

    setUniformArray2fv(name, array) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform2fv(location, array);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform2fv: ${name}`, array);
        }
    }

    setUniformArray3fv(name, array) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform3fv(location, array);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform3fv: ${name}`, array);
        }
    }

    setUniformArray4fv(name, array) {
        const location = this._getUniformLocation(name);
        if (location !== null) {
            this.gl.uniform4fv(location, array);
            if (this.debug) console.debug(`[ShaderBuilder] Set uniform4fv: ${name}`, array);
        }
    }

    _getUniformLocation(name) {
        if (!this.program) throw new Error("Shader program not built.");

        if (!(name in this._uniformCache)) {
            this._uniformCache[name] = this.gl.getUniformLocation(this.program, name);
        }

        const location = this._uniformCache[name];
        if (location === null && this.debug) {
            console.warn(`[ShaderBuilder] Uniform '${name}' not found or not active.`);
        }

        return location;
    }

    // Texture management
    createTexture(name, width, height, format = null, type = null, data = null) {
        const gl = this.gl;
        format = format || gl.RGBA;
        type = type || gl.UNSIGNED_BYTE;

        const texture = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, texture);
        gl.texImage2D(gl.TEXTURE_2D, 0, format, width, height, 0, format, type, data);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

        this._textures.set(name, { texture, width, height, format, type });
        if (this.debug) console.debug(`[ShaderBuilder] Texture created: ${name} (${width}x${height})`);

        return texture;
    }

    bindTexture(name, uniformName, unit = null) {
        const gl = this.gl;
        const textureData = this._textures.get(name);
        if (!textureData) throw new Error(`Texture '${name}' not found.`);

        if (unit === null) {
            unit = this._activeTextureUnit++;
        }

        gl.activeTexture(gl.TEXTURE0 + unit);
        gl.bindTexture(gl.TEXTURE_2D, textureData.texture);

        if (uniformName) {
            this.setUniform1i(uniformName, unit);
        }

        if (this.debug) console.debug(`[ShaderBuilder] Texture bound: ${name} to unit ${unit}`);
        return unit;
    }

    // Framebuffer management for multi-pass rendering
    createFramebuffer(name, width, height, attachments = ['color']) {
        const gl = this.gl;
        const framebuffer = gl.createFramebuffer();
        gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);

        const fbData = {
            framebuffer,
            width,
            height,
            attachments: {},
            colorAttachments: []
        };

        attachments.forEach((attachment, index) => {
            if (attachment === 'color') {
                const colorTexture = gl.createTexture();
                gl.bindTexture(gl.TEXTURE_2D, colorTexture);
                gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, width, height, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
                gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
                gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
                gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
                gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

                gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0 + index, gl.TEXTURE_2D, colorTexture, 0);
                fbData.attachments[`color${index}`] = colorTexture;
                fbData.colorAttachments.push(gl.COLOR_ATTACHMENT0 + index);
            } else if (attachment === 'depth') {
                const depthTexture = gl.createTexture();
                gl.bindTexture(gl.TEXTURE_2D, depthTexture);
                gl.texImage2D(gl.TEXTURE_2D, 0, gl.DEPTH_COMPONENT24, width, height, 0, gl.DEPTH_COMPONENT, gl.UNSIGNED_INT, null);
                gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
                gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
                gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
                gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

                gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.TEXTURE_2D, depthTexture, 0);
                fbData.attachments.depth = depthTexture;
            }
        });

        // Set draw buffers for multiple color attachments
        if (fbData.colorAttachments.length > 0) {
            gl.drawBuffers(fbData.colorAttachments);
        }

        const status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
        if (status !== gl.FRAMEBUFFER_COMPLETE) {
            throw new Error(`Framebuffer '${name}' is not complete: ${status}`);
        }

        this._framebuffers.set(name, fbData);
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);

        if (this.debug) console.debug(`[ShaderBuilder] Framebuffer created: ${name} (${width}x${height})`);
        return framebuffer;
    }

    bindFramebuffer(name) {
        const gl = this.gl;
        if (name === null) {
            gl.bindFramebuffer(gl.FRAMEBUFFER, null);
            if (this.debug) console.debug('[ShaderBuilder] Bound to default framebuffer');
            return;
        }

        const fbData = this._framebuffers.get(name);
        if (!fbData) throw new Error(`Framebuffer '${name}' not found.`);

        gl.bindFramebuffer(gl.FRAMEBUFFER, fbData.framebuffer);
        gl.viewport(0, 0, fbData.width, fbData.height);

        if (this.debug) console.debug(`[ShaderBuilder] Framebuffer bound: ${name}`);
    }

    getFramebufferTexture(framebufferName, attachmentName) {
        const fbData = this._framebuffers.get(framebufferName);
        if (!fbData) throw new Error(`Framebuffer '${framebufferName}' not found.`);

        const texture = fbData.attachments[attachmentName];
        if (!texture) throw new Error(`Attachment '${attachmentName}' not found in framebuffer '${framebufferName}'.`);

        return texture;
    }

    // UBO methods (unchanged)
    createUBO(name, bindingIndex, dataSize, usage = this.gl.DYNAMIC_DRAW) {
        const gl = this.gl;
        const blockIndex = gl.getUniformBlockIndex(this.program, name);
        if (blockIndex === gl.INVALID_INDEX || blockIndex === 4294967295) throw new Error(`Uniform block '${name}' not found.`);
        gl.uniformBlockBinding(this.program, blockIndex, bindingIndex);
        const buffer = gl.createBuffer();
        gl.bindBuffer(gl.UNIFORM_BUFFER, buffer);
        gl.bufferData(gl.UNIFORM_BUFFER, dataSize, usage);
        gl.bindBufferBase(gl.UNIFORM_BUFFER, bindingIndex, buffer);
        this._ubos.set(name, { buffer, bindingIndex, size: dataSize });
        if (this.debug) console.debug(`[ShaderBuilder] UBO created: ${name} bound to ${bindingIndex}`);
    }

    updateUBO(name, data, offset = 0) {
        const gl = this.gl;
        const ubo = this._ubos.get(name);
        if (!ubo) throw new Error(`UBO '${name}' not found.`);
        gl.bindBuffer(gl.UNIFORM_BUFFER, ubo.buffer);
        gl.bufferSubData(gl.UNIFORM_BUFFER, offset, data);
        if (this.debug) console.debug(`[ShaderBuilder] UBO updated: ${name} offset ${offset}`, data);
    }

    bindUBO(name) {
        const gl = this.gl;
        const ubo = this._ubos.get(name);
        if (!ubo) throw new Error(`UBO '${name}' not found.`);
        gl.bindBufferBase(gl.UNIFORM_BUFFER, ubo.bindingIndex, ubo.buffer);
    }

    // VAO methods (unchanged)
    createVAO(name) {
        const vao = this.gl.createVertexArray();
        this._vaos.set(name, vao);
        if (this.debug) console.debug(`[ShaderBuilder] VAO created: ${name}`);
        return vao;
    }

    bindVAO(name) {
        const vao = this._vaos.get(name);
        if (!vao) throw new Error(`VAO '${name}' not found.`);
        this.gl.bindVertexArray(vao);
        if (this.debug) console.debug(`[ShaderBuilder] VAO bound: ${name}`);
    }

    // VBO methods (unchanged)
    createVBO(name, data, usage = this.gl.STATIC_DRAW) {
        const buffer = this.gl.createBuffer();
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer);
        this.gl.bufferData(this.gl.ARRAY_BUFFER, data, usage);
        this._vbos.set(name, buffer);
        if (this.debug) console.debug(`[ShaderBuilder] VBO created: ${name}`);
    }

    setAttribute(name, size, type, normalized = false, stride = 0, offset = 0) {
        const gl = this.gl;
        const location = gl.getAttribLocation(this.program, name);
        if (location === -1) throw new Error(`Attribute '${name}' not found.`);
        gl.enableVertexAttribArray(location);
        gl.vertexAttribPointer(location, size, type, normalized, stride, offset);
        this._attributeCache[name] = location;
        if (this.debug) console.debug(`[ShaderBuilder] Attribute set: ${name}`, { size, type, normalized, stride, offset });
    }

    // Utility methods
    resetTextureUnits() {
        this._activeTextureUnit = 0;
    }

    getTextureInfo(name) {
        return this._textures.get(name);
    }

    getFramebufferInfo(name) {
        return this._framebuffers.get(name);
    }

    // Resource cleanup
    dispose() {
        const gl = this.gl;

        if (this.program) {
            gl.deleteProgram(this.program);
            this.program = null;
        }

        this._uniformCache = {};
        this._attributeCache = {};

        for (const ubo of this._ubos.values()) gl.deleteBuffer(ubo.buffer);
        this._ubos.clear();

        for (const vbo of this._vbos.values()) gl.deleteBuffer(vbo);
        this._vbos.clear();

        for (const vao of this._vaos.values()) gl.deleteVertexArray(vao);
        this._vaos.clear();

        for (const textureData of this._textures.values()) gl.deleteTexture(textureData.texture);
        this._textures.clear();

        for (const fbData of this._framebuffers.values()) {
            gl.deleteFramebuffer(fbData.framebuffer);
            for (const attachment of Object.values(fbData.attachments)) {
                gl.deleteTexture(attachment);
            }
        }
        this._framebuffers.clear();

        if (this.debug) console.debug('[ShaderBuilder] All resources disposed.');
    }

    use() {
        if (!this.program) throw new Error("Shader program not built.");
        this.gl.useProgram(this.program);
        this.resetTextureUnits();
    }
};

/**
 * RenderPipeline: Manages sequential rendering operations with different shaders and targets
 */
class RenderPipeline {
    constructor(gl, debug = false) {
        this.gl = gl;
        this.debug = debug;
        this.steps = [];
        this.sharedUniforms = new Map();
    }

    addStep(name, shader, options = {}) {
        const step = {
            name,
            shader,
            renderTo: options.renderTo || null,
            clearColor: options.clearColor || [0, 0, 0, 1],
            clearDepth: options.clearDepth || 1.0,
            enabled: options.enabled !== false,
            uniforms: new Map(),
            beforeRender: options.beforeRender || null,
            afterRender: options.afterRender || null
        };

        this.steps.push(step);
        if (this.debug) console.debug(`[RenderPipeline] Step added: ${name}`);
        return step;
    }

    enableStep(name, enabled = true) {
        const step = this.steps.find(s => s.name === name);
        if (!step) throw new Error(`Step '${name}' not found.`);
        step.enabled = enabled;
    }

    setStepUniform(stepName, uniformName, value) {
        const step = this.steps.find(s => s.name === stepName);
        if (!step) throw new Error(`Step '${stepName}' not found.`);
        step.uniforms.set(uniformName, value);
    }

    setSharedUniform(name, value) {
        this.sharedUniforms.set(name, value);
    }

    setStepCallbacks(stepName, beforeRender, afterRender) {
        const step = this.steps.find(s => s.name === stepName);
        if (!step) throw new Error(`Step '${stepName}' not found.`);
        step.beforeRender = beforeRender;
        step.afterRender = afterRender;
    }

    render(drawCallback = null) {
        const gl = this.gl;

        for (const step of this.steps) {
            if (!step.enabled) continue;

            // Set render target
            if (step.renderTo) {
                step.shader.bindFramebuffer(step.renderTo);
            } else {
                step.shader.bindFramebuffer(null);
                gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
            }

            // Clear buffers
            gl.clearColor(...step.clearColor);
            gl.clearDepth(step.clearDepth);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            // Activate shader
            step.shader.use();

            // Apply shared uniforms
            for (const [name, value] of this.sharedUniforms) {
                this._applyUniform(step.shader, name, value);
            }

            // Apply step-specific uniforms
            for (const [name, value] of step.uniforms) {
                this._applyUniform(step.shader, name, value);
            }

            // Before render callback
            if (step.beforeRender) {
                step.beforeRender(step, step.shader);
            }

            // Execute drawing
            if (drawCallback) {
                drawCallback(step, step.shader);
            }

            // After render callback
            if (step.afterRender) {
                step.afterRender(step, step.shader);
            }

            if (this.debug) console.debug(`[RenderPipeline] Step rendered: ${step.name}`);
        }
    }

    _applyUniform(shader, name, value) {
        if (Array.isArray(value)) {
            switch (value.length) {
                case 1: shader.setUniform1f(name, value[0]); break;
                case 2: shader.setUniform2f(name, value[0], value[1]); break;
                case 3: shader.setUniform3f(name, value[0], value[1], value[2]); break;
                case 4: shader.setUniform4f(name, value[0], value[1], value[2], value[3]); break;
                case 9: shader.setUniformMatrix3fv(name, value); break;
                case 16: shader.setUniformMatrix4fv(name, value); break;
                default: shader.setUniformArray1fv(name, value); break;
            }
        } else if (typeof value === 'number') {
            if (Number.isInteger(value)) {
                shader.setUniform1i(name, value);
            } else {
                shader.setUniform1f(name, value);
            }
        } else if (typeof value === 'object' && value.type === 'texture') {
            shader.bindTexture(value.name, name, value.unit);
        }
    }

    getStep(name) {
        return this.steps.find(s => s.name === name);
    }

    removeStep(name) {
        const index = this.steps.findIndex(s => s.name === name);
        if (index !== -1) {
            this.steps.splice(index, 1);
            if (this.debug) console.debug(`[RenderPipeline] Step removed: ${name}`);
        }
    }

    dispose() {
        this.steps = [];
        this.sharedUniforms.clear();
        if (this.debug) console.debug('[RenderPipeline] Pipeline disposed.');
    }
}

export { RenderPipeline, ShaderBuilder };
