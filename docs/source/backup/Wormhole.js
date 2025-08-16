import easingUtils from "https://esm.sh/easing-utils";

class BlackHole extends HTMLElement {
  /**
   * Init
   */
  connectedCallback() {
    // Elements
    this.canvas = this.querySelector(".js-canvas");
    this.ctx = this.canvas.getContext("2d");

    // Init
    this.setSizes();

    this.bindEvents();

    // RAF
    requestAnimationFrame(this.tick.bind(this));
  }

  /**
   * Bind events
   */
  bindEvents() {
    window.addEventListener("resize", this.onResize.bind(this));
  }

  /**
   * Resize handler
   */
  onResize() {
    this.setSizes()
  }

  /**
   * Set sizes
   */
  setSizes () {
    this.setCanvasSize()
    this.setGraphics()
  }

  /**
   * Set canvas size
   */
  setCanvasSize () {
    const rect = this.getBoundingClientRect()

    this.render = {
      width: rect.width,
      hWidth: rect.width * 0.5,
      height: rect.height,
      hHeight: rect.height * 0.5,
      dpi: window.devicePixelRatio
    }

    this.canvas.width = this.render.width * this.render.dpi
    this.canvas.height = this.render.height * this.render.dpi
  }

  /**
   * Set graphics
   */
  setGraphics () {
    this.setDiscs()
    this.setDots()
  }

  /**
   * Set discs
   */
  setDiscs () {
    this.discs = []

    this.startDisc = {
      x: this.render.width * 0.5,
      y: this.render.height * 0,
      w: this.render.width * 1,
      h: this.render.height * 1
    }

    const totalDiscs = 150

    for (let i = 0; i < totalDiscs; i++) {
      const p = i / totalDiscs

      const disc = this.tweenDisc({
        p
      })

      this.discs.push(disc)
    }
  }

  /**
   * Set dots
   */
  setDots () {
    this.dots = []

    const totalDots = 20000

    for (let i = 0; i < totalDots; i++) {
      const disc = this.discs[Math.floor(this.discs.length * Math.random())]
      const dot = {
        d: disc,
        a: 0,
        c: `rgb(${Math.random() * 0}, ${150 + Math.random() * 50}, ${150 + Math.random() * 105})`,
        p: Math.random(),
        o: Math.random()
      }

      this.dots.push(dot)
    }
  }


  /**
   * Tween disc
   */
  tweenDisc (disc) {
    const { startDisc } = this

    const scaleX = this.tweenValue(1, 0, disc.p, 'outCubic')
    const scaleY = this.tweenValue(1, 0, disc.p, 'outExpo')

    disc.sx = scaleX
    disc.sy = scaleY

    disc.w = startDisc.w * scaleX
    disc.h = startDisc.h * scaleY

    disc.x = startDisc.x
    disc.y = startDisc.y + disc.p * startDisc.h * 1

    return disc
  }

  /**
   * Tween value
   */
  tweenValue (start, end, p, ease = false) {
    const delta = end - start

    const easeFn =
      easingUtils[
        ease ? "ease" + ease.charAt(0).toUpperCase() + ease.slice(1) : "linear"
      ];

    return start + delta * easeFn(p);
  }

  /**
   * Draw discs
   */
  drawDiscs () {
    const { ctx } = this

    ctx.strokeStyle = '#0329'
    ctx.lineWidth = 1

    // Discs
    this.discs.forEach((disc, i) => {
      const p = disc.sx * disc.sy

      ctx.beginPath()

      ctx.globalAlpha = disc.a

      ctx.ellipse(
        disc.x,
        disc.y + disc.h,
        disc.w,
        disc.h,
        0,
        0,
        Math.PI * 2
      )
      ctx.stroke()

      ctx.closePath()
    })
  }

  /**
   * Draw dots
   */
  drawDots () {
    const { ctx } = this


    this.dots.forEach((dot) => {
      const { d, a, p, c, o } = dot

      const _p = d.sx * d.sy
      ctx.fillStyle = c

      const newA = a + (Math.PI * 2 * p)
      const x = d.x + Math.cos(newA) * d.w
      const y = d.y + Math.sin(newA) * d.h

      ctx.globalAlpha = d.a * o

      ctx.beginPath()
      ctx.arc(x, y + d.h, 1 + _p * 0.5, 0, Math.PI * 2)
      ctx.fill()
      ctx.closePath()
    })
  }

  /**
   * Move discs
   */
  moveDiscs () {
    this.discs.forEach((disc) => {
      disc.p = (disc.p + 0.0003) % 1

      this.tweenDisc(disc)

      const p = disc.sx * disc.sy

      let a = 1
      if (p < 0.01) {
        a = Math.pow(Math.min(p / 0.01, 1), 3)
      } else if (p > 0.2) {
        a = 1 - Math.min((p - 0.2) / 0.8, 1)
      }

      disc.a = a
    })
  }

  /**
   * Move dots
   */
  moveDots () {
    this.dots.forEach((dot) => {
      const v = this.tweenValue(0, 0.001, 1 - dot.d.sx * dot.d.sy, 'inExpo')
      dot.p = (dot.p + v) % 1
    })
  }

  /**
   * Tick
   */
  tick(time) {
    const { ctx } = this;

    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

    ctx.save();
    ctx.scale(this.render.dpi, this.render.dpi);

    // Move
    this.moveDiscs()
    this.moveDots()

    // Draw
    this.drawDiscs()
    this.drawDots()

    ctx.restore();

    requestAnimationFrame(this.tick.bind(this));
  }
}

customElements.define("black-hole", BlackHole);
