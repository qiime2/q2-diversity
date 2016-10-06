/* global d */
/* global document */
/* global XMLSerializer */

import { select } from 'd3';

import setupData from './data';
import { render, kwStats, warnings } from './render';


export function addCategoryPicker(row, categories, selectedCategory) {
  const grp = row.append('div').attr('class', 'col-lg-2 form-group categoryPicker');
  grp.append('label').text('Category');
  grp.append('select')
    .attr('class', 'form-control')
    .on('change', function changeCategory() {
      const data = d[this.selectedIndex];
      const preppedData = setupData(data);
      const svg = select('svg');
      render(svg, preppedData);
      const body = select('body .container-fluid');
      kwStats(body, data);
      warnings(body, data);
    })
    .selectAll('option')
    .data(categories)
    .enter()
      .append('option')
      .attr('value', d => d)
      .text(d => d)
      .property('selected', d => (d === selectedCategory));
  return grp;
}

export function addDownloadLinks(sel, svg) {
  const grp = sel.append('div').attr('class', 'col-lg-2 form-group');
  grp.append('label').html('&nbsp;');
  grp.append('button')
    .text('Download SVG')
    .attr('class', 'btn btn-default form-control')
    .on('click', () => {
      const serializer = new XMLSerializer();
      let src = serializer.serializeToString(svg.node());
      src = `<?xml version="1.0" standalone="no"?>\r\n${src}`;
      const url = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(src)}`;
      const link = document.createElement('a');
      link.setAttribute('href', url);
      link.setAttribute('download', 'alpha-compare.svg');
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    });
}
