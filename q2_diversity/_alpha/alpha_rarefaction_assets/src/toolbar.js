/* global d */
/* global document */
/* global XMLSerializer */

import { select } from 'd3';

// import setupData from './data';
// import { render, stats, warnings } from './render';

export function addMetricPicker(row, metrics, selectedMetric) {
  const grp = row.append('div').attr('class', 'col-lg-2 form-group metricPicker');
  grp.append('label').text('Metric');
  grp.append('select')
    .attr('class', 'form-control')
    .selectAll('option')
    .data(metrics)
    .enter()
      .append('option')
      .attr('value', d => d)
      .text(d => d)
      .property('selected', d => (d === selectedMetric));
  return grp;
}

export function addDownloadLinks(sel) {
  const grp = sel.append('div')
    .attr('class', 'col-lg-2 form-group')
    .attr('id', 'downloadLink');
  grp.append('label').html('&nbsp;');
  grp.append('button')
    .text('Download CSV')
    .attr('class', 'btn btn-default form-control')
    .on('click', () => {
      const link = document.createElement('a');
      const selectedMetric = select('form-control').node().value;
      const url = `metric-${selectedMetric}.csv`;
      link.setAttribute('href', url);
      link.setAttribute('download', url);
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    });
}
