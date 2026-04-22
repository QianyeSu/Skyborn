document.addEventListener("DOMContentLoaded", () => {
  const isTutorialIndexLink = (anchor) => {
    const href = anchor.getAttribute("href");
    if (!href) {
      return false;
    }

    try {
      const resolved = new URL(href, window.location.href);
      return resolved.pathname.endsWith("/notebooks/index.html");
    } catch (_error) {
      return false;
    }
  };

  const collapseTutorialSidebarBranch = () => {
    const sidebarLinks = document.querySelectorAll(
      ".bd-sidebar-primary .bd-sidenav a.reference.internal"
    );

    sidebarLinks.forEach((anchor) => {
      if (!isTutorialIndexLink(anchor)) {
        return;
      }

      const item = anchor.closest("li");
      if (!item) {
        return;
      }

      item.classList.add("skyborn-tutorials-root");
      item.classList.remove("has-children");

      const details = item.querySelector(":scope > details");
      if (details) {
        details.removeAttribute("open");
      }
    });
  };

  collapseTutorialSidebarBranch();
  window.setTimeout(collapseTutorialSidebarBranch, 150);
});
