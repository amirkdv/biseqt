$(document).ready(function() {
  $('img[alt]').each(function() {
    if($(this).attr('alt') == 'lightbox'){
      a_tag = $(this).parent()
      a_tag.attr('rel', 'lightbox')

      caption = a_tag.next()
      if(caption.is('p.caption')) {
        a_tag.attr('data-title', caption.text())
      }
    }
  });
});
